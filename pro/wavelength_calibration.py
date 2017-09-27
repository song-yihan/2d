import numpy as np
from scipy import optimize
import logging
import extract_1d
import os
import astropy.io.fits as pyfits

class wavelength_calibration():
    def __init__(self, arc_spec=None):
        self.wavelength_points = {}
        self.wavelength = {}
        self.arc = {}
        self.peaks = {}
        self.line_table = []
        self.arc_spec = arc_spec
        if arc_spec is None or not isinstance(arc_spec, extract_1d.extract_1d):
            #arc_spec = extract_1d.extract_1d()
            return
        else:
            for spec in arc_spec:
                fiberid = spec.fiberid
                flux = spec.flux
                self.arc[fiberid] = flux
            #self.arc_spec = arc_spec
        pass
    def clean(self):
        self.wavelength_points = {}
        self.wavelength = {}
        self.arc = {}
        self.peaks = {}
        self.line_table = []
    def load_arc_lines(self, ff):
        lines = {}
        for line in open(ff):
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            line = line.replace(',',' ')
            line = line.replace('\t',' ')
            while '  ' in line:
                line = line.replace('  ',' ')
            items = line.split(' ')
            if len(items) !=3:
                continue
            lines[int(items[0])] = {'wavelength':float(items[1]),'name':items[2]}
        self.line_table = lines
    def set_arc(self, arc):
        self.arc = {}
        for i in range(len(arc)):
            self.arc[i] = arc[i]
    def make_wavelength(self, fiberid = 'all', PN=5):
        if fiberid == 'all':
            fiberids = self.wavelength_points.keys()
        else:
            fiberids = [fiberid]
        if len(fiberids) == 0:
            return
        fiberids.sort()
        n = self.wavelength_points[fiberids[0]]
        x = np.array(range(len(self.arc[fiberids[0]])))
        for fiberid in fiberids:
            wp = self.wavelength_points[fiberid]
            if len(wp['center']) +1 <= PN:
                logging.debug('Number of the detected lines is too small to fit. Fiberid:%d'%(fiberid))
                continue
            coeff = np.polyfit(wp['center'], wp['wavelength'], PN)
            wl = np.polyval(coeff, x)
            self.wavelength[fiberid] = {'coeff':coeff, 'wavelength':wl}
    def add_wavelength_points(self, fiberid, peak_ind, line_table_ind):
        if not self.wavelength_points.has_key(fiberid):
            self.wavelength_points[fiberid] = {'center':[], 'wavelength':[],'index':[]}
        if line_table_ind >= len(self.line_table) or line_table_ind < 0:
            return
        if not self.peaks.has_key(fiberid):
            return
        peaks = self.peaks[fiberid]
        #print peaks
        peak_x = peaks['wavelength']
        peak_index = peaks['index']
        if peak_ind <0 or peak_ind >= len(peak_x):
            return

        center = peak_x[peak_ind]
        index = peak_index[peak_ind]
        wavelength = self.line_table[line_table_ind]['wavelength']
        wp = self.wavelength_points[fiberid]
        n = len(wp['center'])
        i = 0
        while i<n and wp['center'][i] < center:
            i+=1
        if i == n:
            wp['center'].append(center)
            wp['wavelength'].append(wavelength)
            wp['index'].append(index)
            return
        wp['center'] = wp['center'][:i]+[center]+wp['center'][i:]
        wp['wavelength'] = wp['wavelength'][:i]+[wavelength]+wp['wavelength'][i:]
        wp['index'] = wp['index'][:i]+[index]+wp['index'][i:]

    def clean_wavelenth_points(self, fiberid):
        if not self.wavelength_points.has_key(fiberid):
            return
        wp = self.wavelength_points[fiberid]
        while len(wp['center']) != 0:
            del wp['center'][0]
            del wp['wavelength'][0]
            del wp['index'][0]

    def block_peak(self, fiberid, wavelength_range):
        wave_start, wave_end = wavelength_range
        if not self.peaks.has_key(fiberid) or not self.wavelength.has_key(fiberid):
            return
        coeff = self.wavelength[fiberid]['coeff']
        peaks = self.peaks[fiberid]
        peak_wavelength = peaks['wavelength']
        i = 0
        while i < len(peak_wavelength):
            #print i, len(peak_wavelength)
            w = peak_wavelength[i]
            angstrom = np.polyval(coeff, w)
            if angstrom < wave_start or angstrom > wave_end:
                i += 1
                continue
            del peaks['index'][i]
            del peaks['flux'][i]
            del peaks['wavelength'][i]

    def match_lines_sc_blue(self, fiberid):
        if not self.peaks.has_key(fiberid):
            return
        peaks = self.peaks[fiberid]
        #print peaks
        x = peaks['wavelength']
        ind = 0
        while ind < len(x) and x[ind]<2000:
            ind+=1
        #print ind
        if ind == len(x):
            return

        self.clean_wavelenth_points(fiberid)
        self.add_wavelength_points(fiberid,ind-1,6)
        self.add_wavelength_points(fiberid,ind,7)
        self.add_wavelength_points(fiberid,ind+1,8)
        self.add_wavelength_points(fiberid,ind+2,9)
        self.add_wavelength_points(fiberid,ind+4,10)
        self.add_wavelength_points(fiberid,ind+5,11)
        self.add_wavelength_points(fiberid,ind+6,12)
        self.make_wavelength(fiberid, PN=2)
        self.block_peak(fiberid, (5080, 5087.74))
        ind = 0
        while ind < len(x) and x[ind]<2000:
            ind+=1
        #print ind
        if ind == len(x):
            return

        self.double_line_fit(fiberid, (ind-16,ind-15))
        self.add_wavelength_points(fiberid,ind-1,6)
        self.add_wavelength_points(fiberid,ind-8,5)
        self.add_wavelength_points(fiberid,ind-9,4)
        self.add_wavelength_points(fiberid,ind-10,3)
        self.add_wavelength_points(fiberid,ind-15,2)
        self.add_wavelength_points(fiberid,ind-16,1)
        self.add_wavelength_points(fiberid,ind,7)
        self.add_wavelength_points(fiberid,ind+1,8)
        self.add_wavelength_points(fiberid,ind+2,9)
        self.add_wavelength_points(fiberid,ind+4,10)
        self.add_wavelength_points(fiberid,ind+5,11)
        self.add_wavelength_points(fiberid,ind+6,12)
    def double_line_fit(self, fiberid, pair_index):
        peaks = self.peaks[fiberid]
        peak_x = peaks['wavelength']
        peak_index = peaks['index']
        p1, p2 = pair_index
        p1_index, p2_index = peak_index[p1], peak_index[p2]
        y = self.arc[fiberid]
        x = np.array(range(len(y)))
        half_width = 10
        vec_x = x[p1_index-half_width:p2_index+half_width+1]
        vec_y = y[p1_index-half_width:p2_index+half_width+1]
        vec1 = y[p1_index-half_width:p1_index+half_width+1]
        vec2 = y[p2_index-half_width:p2_index+half_width+1]
        mag1 = np.max(vec1)-np.min(vec1)
        mag2 = np.max(vec2)-np.min(vec2)
        p0 = [mag1,peak_x[p1],3.5,mag2,peak_x[p2],3.5,np.min(vec_y)]
        gauss_func = lambda p, x: p[0]*np.exp(-0.5*(x-p[1])**2/(p[2]**2))+p[3]*np.exp(-0.5*(x-p[4])**2/(p[5]**2))+p[6]
        err_func = lambda p, x,y: gauss_func(p,x)- y
        p, success = optimize.leastsq( err_func, p0, args=(vec_x, vec_y))
        if not success in [1,2,3,4]:
            return

        peak_x[p1] = p[1]
        peak_x[p2] = p[4]
        #plt.plot(vec_x, vec_y,'r')
        pass
    def find_peak(self,fiberid, half_width = 10):
        arc = self.arc[fiberid]
        
        y = arc
        x = np.array(range(len(y)))
        gap = half_width
        ind1 = y[1:] >= y[:-1]
        ind2 = y[:-1] >= y[1:]
        ind3 = ind1[:-1] & ind2[1:]
        xx = x[1:-1][ind3]
        yy = y[1:-1][ind3]
        #print len(x),len(ind3),len(d),len(y)
        xx,yy  = list(xx), list(yy)
        #print xx
        i = 0
        while i<len(xx)-1:
            if xx[i+1] - xx[i] > gap:
                i+=1
                continue
            del_index = -1
            if yy[i] >= yy[i+1]:
                del_index = i+1
            else:
                del_index = i
            del xx[del_index]
            del yy[del_index]
        while xx[0]<half_width:
            del xx[0]
            del yy[0]
        n = len(x)
        while xx[-1]>=n-half_width-1:
            del xx[-1]
            del yy[-1]


        ew_list = []

        for i in range(len(xx)):
            ind = xx[i]
            vec_x = x[ind-half_width:ind+half_width+1]
            vec_y = y[ind-half_width:ind+half_width+1]
            ew = self.compute_ew(vec_x, vec_y)
            ew_list.append(ew)
        largest_ew_value = np.max(ew_list)
        valid_min_ew_value = largest_ew_value*0.1
        i = 0
        while i<len(xx):
            if ew_list[i] < valid_min_ew_value:
                del ew_list[i]
                del xx[i]
                del yy[i]
            else:
                i+=1
        i = 0
        out_center = []
        while i<len(xx):
            ind = xx[i]
            vec_x = x[ind-half_width:ind+half_width+1]
            vec_y = y[ind-half_width:ind+half_width+1]
            result = self.gaussian_fit(vec_x, vec_y)
            if result is None or result[0] < 0 or np.abs(result[2]<1.5) or np.abs(result[2]>8):
                del xx[i]
                del yy[i]
                del ew_list[i]
            else:
                out_center.append(result[1])
                i+=1
        out_x = xx
        out_y = yy
        
        #out_center = np.array(out_center)
        self.peaks[fiberid] = {'index':out_x,'flux': out_y,'wavelength':out_center}
        return np.array(out_x), np.array(out_y)
    def compute_ew(self, xx,yy):
        x, y = np.array(xx), np.array(yy)
        n = len(x)
        hn = int(n/2)-1
        #left_ind = y[:hn].argmin()
        #right_ind = y[-hn:].argmin()
        #right_ind += n-hn
        #x0,x1 = x[left_ind],x[right_ind]
        #y0,y1 = y[left_ind],y[right_ind]
        x0, x1 = x[0], x[-1]
        y0, y1 = y[0], y[-1]
        baseline_coeff = np.polyfit([x0,x1],[y0,y1],1)
        baseline = np.polyval(baseline_coeff,x)
        y = y - baseline
        value = (y[:-1]+y[1:])*(x[1:]-x[:-1])*0.5
        value = np.sum(value)
        
        return value
    def gaussian_fit(self, xx, yy):
        center = np.sum(xx * yy) / np.sum(yy)
        mag = yy.max() - yy.min()
        baseline = yy.min()
        p0 = [mag, center, 2.5, baseline]
        gauss_func = lambda p, x: p[0]*np.exp(-0.5*(x-p[1])**2/(p[2]**2))+p[3]
        err_func = lambda p, x,y: gauss_func(p,x)- y
        p, success = optimize.leastsq( err_func, p0, args=(xx, yy))
        if success not in [1,2,3,4]:
            peak = center
            return None
        return p
        pass
    def load(self, filename):
        if not os.path.exists(filename):
            logging.info('\tWavelength load ERROR. File not exists.\n%s'%(filename))

            return
    def save(self, filename):
        pass
    def __len__(self):
        return len(self.wavelength.keys())
    def __getitem__(self,k):
        return self.wavelength[k]['wavelength']
    def __iter__(self):
        i = 0 
        while True:
            try:
                yield self[i]
                i += 1
            except IndexError:
                break
    def has_key(self,key):
        return self.wavelength.has_key(key) 
    def calibration_sc_blue(self):
        self.load_arc_lines('linetable_sc.txt')
        fiberids = self.arc.keys()
        fiberids.sort()
        
        for fiberid in fiberids:
            logging.debug('Wavelength Calibration Fiberid: %d'%(fiberid))
            if self.arc_spec is not None and self.arc_spec[fiberid].tag !=0 :
                # the used trace line is not from trace. That is from inherit
                # the fiber is not good
                continue
            self.find_peak(fiberid)        
            logging.debug('Number of peaks:%d'%(len(self.peaks[fiberid]['flux'])))
            if len(self.peaks[fiberid]['flux'])>27:
                continue
            self.match_lines_sc_blue(fiberid)
        self.make_wavelength(fiberid = 'all')

    

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import astropy.io.fits as pyfits
    import trace2d
    import image2d
    import extract_1d

    #hdulist = pyfits.open('test/02b_arc.fits')
    #data = hdulist[0].data

    trace_data = trace2d.trace()
    trace_data.load('output_spec/trace_footprint.fits')
    
    arc_image = image2d.image2d('../data/20170610/arc/ra-02b-20170610184326-2749-83397283.fit.gz')
    arc_image.subtract_overscan()
    arc_data = extract_1d.extract_1d(arc_image)
    arc_data.set_trace(trace_data)
    arc_data.extract_all()
    #arc_data.writeto('test_arc.fits')

    wc = wavelength_calibration(arc_data)
    wc.calibration_sc_blue()

    obj = image2d.image2d('../data/20170610/BD47293601/ro-02b-20170611024033-16-83397760.fit.gz')
    obj.subtract_overscan()

    data = extract_1d.extract_1d(obj)
    data.set_trace(trace_data)
    data.set_wavelength(wc)
    data.extract_all()
    data.writeto('test.fits')
