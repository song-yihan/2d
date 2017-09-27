import image2d
import trace2d
import numpy as np
import os
import astropy.io.fits as pyfits

class extract_1d():
    def __init__(self, image=None):
        self.image = image
        if not isinstance(image, image2d.image2d):
            self.image = image2d.image2d()
        self.spectra = {}
        self.kernel = {}
        self.trace_lines = []
        self.wavelength_data = {}
    '''
    def trace(self, calibration_data = None):   
        image = self.image
        if image.is_valid():
            
            trace_lines = trace2d.trace(image)
            #print len(lines)
            trace_lines.calibration(calibration_data)
            self.calibration_data = calibration_data
            #print len(lines)
        else:
            trace_lines = trace2d.trace()
            self.calibration_data = None

        self.trace_lines = trace_lines
    '''
    def set_trace(self,t):

        self.trace_lines = t
    def set_wavelength(self, wc):
        self.wavelength_data = wc
    def set_image(self,img):
        self.image = img
    def set_kernel(self,k):
        self.kernel = k
    def extract_all(self):
        image = self.image
        if not image.is_valid():
            return


        for fiberid in range(len(self.trace_lines)):
            
            self.extract_spectra_sum(fiberid)
    def extract_all_kernel(self):
        image = self.image
        if not image.is_valid():
            return


        for fiberid in range(len(self.trace_lines)):
            
            self.extract_spectra_kernel(fiberid)
    def extract_spectra(self, fiberid=0):
        self.extract_spectra_sum(fiberid)


    def extract_spectra_peak(self, fiberid=0):
        image = self.image
        if not image.is_valid():
            return
        line = self.trace_lines[fiberid]

        data = image.get_image()
        fiber_half_width = self.trace_lines.fiber_half_width
        flux = []
        n = len(line['x'])
        for i in range(n):
            index_x, index_y = line['x'][i], line['y'][i]
            flux.append(data[index_x, index_y])

        self.spectra[fiberid] = spectra_1d(flux=np.array(flux), fiberid=fiberid )
    
        
    def extract_spectra_kernel(self, fiberid=0):

        image = self.image
        if not image.is_valid():
            return
        line = self.trace_lines[fiberid]

        data = image.get_image()
        fiber_half_width = self.trace_lines.fiber_half_width
        flux = []
        n = len(line['x'])
        if self.kernel.has_key(fiberid):
            kernel_list = self.kernel[fiberid]
        else:
            kernel_list = []
            for i in range(n):
                index_x, index_y = line['x'][i], line['y'][i]
                y = data[index_x, index_y-fiber_half_width:index_y+fiber_half_width+1]
                x = range(len(y))
                a = self.trace_lines.gaussian_fit_single(x,y)
                a = a[1]
                a = a/np.sum(a)

                kernel_list.append(a)
            self.kernel[fiberid] = kernel_list
        for i in range(n):
            index_x, index_y = line['x'][i], line['y'][i]
            y = data[index_x, index_y-fiber_half_width:index_y+fiber_half_width+1]
            a = kernel_list[i] 
            flux.append(np.sum(y*a))
        self.spectra[fiberid] = spectra_1d(flux=np.array(flux), fiberid=fiberid )
            
    def extract_spectra_sum(self, fiberid=0):
        
        image = self.image
        if not image.is_valid():
            return
        line = self.trace_lines[fiberid]

        data = image.get_image()
        fiber_half_width = self.trace_lines.fiber_half_width
        flux = []
        error = []
        invvar = []
        n = len(line['x'])
        for i in range(n):
            index_x, index_y = line['x'][i], line['y'][i]
            y = data[index_x, index_y-fiber_half_width:index_y+fiber_half_width+1]
            value = np.sum(y)
            flux.append(value)

            y_err = image.error[index_x, index_y-fiber_half_width:index_y+fiber_half_width+1]
            value_err = 0
            for i in range(len(y_err)):
                value_err += y_err[i]**2
            error.append(np.sqrt(value_err))
            invvar.append(1.0 / value_err)

        if self.wavelength_data.has_key(fiberid):
            wave = self.wavelength_data[fiberid]
        else:
            wave = np.array(range(len(flux)))
        tag_str = line['tag']
        if tag_str == 'find':
            tag = 0
        elif tag_str == 'calibration':
            tag = 1
        self.spectra[fiberid] = spectra_1d(wavelength =  wave, flux=np.array(flux), fiberid=fiberid ,tag=tag, error = error, invvar = invvar)
    def writeto_one(self,fiberid=None, filename = None):
        if fiberid is None or filename is None:
            return
    def writeto(self,filename=None):
        self.writeto_all(filename)
    def writeto_all(self,filename=None):
        if filename is None:
            return
        spec = []
        wavelength = []
        tag = []
        invvar = []
        for i in range(len(self)):
            try:
                flux = self[i].flux
                wave = self[i].wavelength
                spec_tag = self[i].tag
                value_invvar = self[i].invvar
            except IndexError:
                continue
            spec.append(flux)
            wavelength.append(wave)
            tag.append(spec_tag)
            invvar.append(value_invvar)
        spec = np.array(spec)
        wavelength = np.array(wavelength)
        tag = np.array(tag)
        invvar = np.array(invvar)

        #print spec.shape
        hdu_spec = pyfits.PrimaryHDU(spec)
        hdu_wavelength = pyfits.ImageHDU(wavelength)
        hdu_tag = pyfits.ImageHDU(tag)
        hdu_invvar = pyfits.ImageHDU(invvar)
        hdulist = pyfits.HDUList([hdu_spec,hdu_wavelength,hdu_invvar,hdu_tag])
        if os.path.exists(filename):
            os.remove(filename)
        hdulist.writeto(filename)

    def __len__(self):
        return len(self.spectra.keys())
    def __getitem__(self,k):
        return self.spectra[k]
    def __iter__(self):
        keys = self.spectra.keys()
        keys.sort()
        i = 0 
        while True:
            try:
                yield self[keys[i]]
                i += 1
            except IndexError:
                break
    def has_key(self,key):
        return self.spectra.has_key(key)

class spectra_1d():
    def __init__(self, wavelength=None, flux=None, fiberid=-1, error=None,invvar = None,tag = 0):
        self.wavelength = np.array(wavelength)
        self.flux = np.array(flux)
        self.fiberid = fiberid
        self.error = np.array(error)
        self.tag = tag
        self.invvar = np.array(invvar)
        pass

if __name__ == '__main__':
    import pickle
    import matplotlib.pyplot as plt
    ff = './tmp/rs-04b-20170414184339-2604-83315203.fit.gz'
    image = image2d.image2d(ff)
    image.subtract_overscan()
    #t = trace.trace(image)
    fiberid = 2
    import time
    trace_file = './tmp/trace.pkl'

    data1d = extract_1d(image)

    t1 = time.time()
    data1d.trace(None)
    #data1d.set_trace(pickle.load(open(trace_file)))
    t1 = time.time()-t1
    
    #trace_f = open(trace_file,'w')
    #pickle.dump(data1d.trace_lines, trace_f)
    #trace_f.close()
    

    t2 = time.time() 
    data1d.extract_spectra_peak(fiberid)
    spectra_peak = data1d[fiberid]
    t2 = time.time()-t2

    t3 = time.time()
    data1d.extract_spectra_kernel(fiberid)
    spectra_kernel = data1d[fiberid]
    t3 = time.time()-t3

    t4 = time.time()
    data1d.extract_spectra_sum(fiberid)
    spectra_sum = data1d[fiberid]
    t4 = time.time()-t4

    print 'trace time spends:',t1
    print 'extrace_peak spends:',t2
    print 'extrace_kernel spends:',t3
    print 'extrace_sum spends:',t4
    flux_peak = spectra_peak.flux
    flux_kernel = spectra_kernel.flux
    flux_sum = spectra_sum.flux

    from scipy import optimize


    p0 = [1]
    err_func = lambda p,x,y:x*p[0] - y
    p,_ = optimize.leastsq(err_func, p0, args=(flux_kernel, flux_peak))
    flux_kernel = flux_kernel * p[0]
    p0 = [1]
    p,_ = optimize.leastsq(err_func, p0, args=(flux_sum, flux_peak))
    flux_sum = flux_sum * p[0]

    import matplotlib.pyplot as plt
    plt.figure(0)
    plt.plot(flux_peak, label='peak')
    plt.plot(flux_kernel, label='kernel')
    plt.plot(flux_sum, label='sum')
    plt.legend(loc='best')
    #plt.figure(1)
    #plt.plot(flux-flux_test)
    #plt.figure(2)
    #plt.plot((flux-flux_test)/flux)
    #plt.figure(3)
    #plt.plot((flux-
    ''' 
    hexx, hexy = [],[]
    plt.figure(4)
    for i in range(len(data1d.kernel)):
        y = data1d.kernel[i]
        x = range(len(y))
        for j in range(len(x)):
            hexx.append(x[j])
            hexy.append(y[j])
        #plt.plot(d,'b')
    plt.hexbin(hexx,hexy,mincnt=1)
    plt.colorbar()
    '''
    plt.figure(5)
    x,y = [],[]
    x1,y1 = [],[]
    kernel_list = data1d.kernel[fiberid]
    for i in range(len(kernel_list)-1):
        x.append(i)
        k1 = kernel_list[i]
        k2 = kernel_list[i+1]
        
        value = np.sqrt(np.sum((k1-k2)**2))
        y.append(value)

        x1.append(i)
        y1.append(np.abs(k2.argmax()-k1.argmax()))
    plt.plot(x,y)
    y1 = np.array(y1)/25
    x1 = np.array(x1)
    plt.plot(x1,y1,'ro-')

    plt.figure(6)
    plt.plot(kernel_list[1546],'r',label='1546')
    plt.plot(kernel_list[1547],'g',label='1547')
    plt.plot(kernel_list[1548],'y',label='1548')
    plt.legend(loc='best') 
    
    plt.show()
