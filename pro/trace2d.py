import numpy as np
from scipy import optimize
import image2d
import matplotlib.pyplot as plt
import pickle
import sys
import os
import astropy.io.fits as pyfits

class trace():
    def __init__(self, image=None, sampling_step = 50):
        self.image = image
        self.fiber = []    
        self.fiber_half_width = 6
        self.sampling_step = sampling_step
        if not isinstance(image, image2d.image2d):
            self.image = image2d.image2d()
            return
        if not image.is_valid():
            return
        self.trace_image()
    def remove_image(self):
        del self.image
        self.image = None
    def calibration(self, cali_array):
        '''
        find the connection between fibers we found and cali_array.
        After this progress, 250 fibers are confirmed.
        :param cali_array: 250 x 4116 float array. Each row stands for one fiber footprint
        :return:
            self.fiber stores 250 fibers' center footprint
        '''
        if isinstance(cali_array, str):
            try:
                hdulist = pyfits.open(cali_array)
                data = hdulist[0].data
            except:
                return
            cali_array = data
        if not isinstance(cali_array, np.ndarray):
            return
        if len(cali_array) == 0:
            return

        n = len(cali_array[0])
        x = np.array(range(n))
        cali_lines = []
        for i in range(len(cali_array)):

            vector = cali_array[i]
            y = np.array([int(k+0.5) for k in vector])
            item = {'x':x[:], 'y':y, 'y_center':vector,'tag':'calibration'}
            cali_lines.append(item)

        compare_index = 2000
        fiber = self.fiber
        fiber_ind, cali_ind = 0,0
        fiber_gap_threshold = 6.0
        while fiber_ind<len(fiber) and cali_ind < len(cali_lines):
            fiber_value = fiber[fiber_ind]['y_center'][compare_index] 
            cali_value  = cali_lines[cali_ind]['y_center'][compare_index]
            d = fiber_value - cali_value
            if np.abs(d) < fiber_gap_threshold:
                # the fiber index matches the calibration line index
                # Each index plus one and go to the next index
                fiber_ind +=1
                cali_ind +=1
                continue
            elif d >=fiber_gap_threshold: 
                ## not matches ##
                ## fiber value is larger than cali value
                ## insert the cali line into the fiber at fiber_index
                fiber = fiber[:fiber_ind]+ [cali_lines[cali_ind]] +fiber[fiber_ind:]
                cali_ind +=1
                fiber_ind +=1

            elif d <=-fiber_gap_threshold:
                ## not matches ##
                ## fiber value is smaller than cali value
                ## fiber index plus one and continue
                fiber_ind +=1
                continue
        if cali_ind < len(cali_lines):
            fiber = fiber + cali_lines[cali_ind:]
        if len(fiber) == 250:
            self.fiber = fiber
        else:
            # error
            self.fiber = cali_lines
        
        return
    def set_trace(self, data):
        '''
        Input fiber centers from outside
        :param data: list
        :return:
            let self.fiber = data
        '''
        self.fiber = data
    def __len__(self):
        return len(self.fiber)
    def __getitem__(self,k):
        return self.fiber[k]
    def __iter__(self):
        i = 0 
        while True:
            try:
                yield self[i]
                i += 1
            except IndexError:
                break

    def find_best_friend(self, vec1, vec2):
        '''
        For each point in vec1, find the closed point in vec2.
        :param vec1: a list of centers
        :param vec2: a list of centers
        :return:
            return the index list. the index is the location of the best frind in the vec2.
        '''
        n = len(vec1)
        m = len(vec2)
        i = 0
        j = 0
        rec = np.zeros(n, dtype=int)
        for i in range(n):
            best_value = np.abs(vec1[i] - vec2[j])
            best_ind = j
            j+=1
            while j<m:
                value = np.abs(vec1[i] - vec2[j])
                if value> best_value:
                    break
                if value<= best_value:
                    best_value = value
                    best_ind = j
                j+=1
            rec[i] = best_ind
            j = best_ind
        return rec
    def peak_cluster_Recursively(self):
        '''
        Called by peak_cluster

        cluster the centers into groups.
        Using the theory: I am my best frind's best frind.
         Because the searching order is from top to bottom and left to the right,
         the list of self.cluster_dict[key]['x'] is from small to large.
        :return:
            self.cluster_dict: Dict, key is the group number,
                                    value:{'x':row number list,'y':peak center list}
        '''
        self.cluster_dict = {}
        for ln in range(len(self.line_number_list)):
            line_key = self.line_number_list[ln]
            for num in range(len(self.cluster_array_dict[line_key])):
                cluster_array_dict = self.cluster_array_dict
                current_line = cluster_array_dict[line_key]
                if current_line[num] != 0:
                    continue
                self.cluster_number +=1
                cn = self.cluster_number
                if self.cluster_dict.has_key(cn):
                    self.cluster_dict[cn]['x'].append(line_key)
                    self.cluster_dict[cn]['y'].append(self.line_peak[line_key][num])
                else:
                    self.cluster_dict[cn] = {'x':[line_key],
                                             'y':[self.line_peak[line_key][num]]
                                             }
                current_line[num] = cn
                '''
                x, y = ln, num
                while x>0:
                    k = self.line_number_list[x]
                    x1, y1 = x-1, self.best_friend2[k][y]
                    k1 = self.line_number_list[x1]
                    if y != self.best_friend1[k1][y1]: # not good friend
                        break
                    self.cluster_array_dict[k1][y1] = cn
                    x, y = x1,y1
                '''
                x, y = ln, num
                while x<len(self.line_number_list)-1:
                    k = self.line_number_list[x]
                    x1, y1 = x+1, self.best_friend1[k][y]
                    k1 = self.line_number_list[x1]
                    if y != self.best_friend2[k1][y1]: #not good friend
                        break
                    self.cluster_array_dict[k1][y1] = cn
                    if self.cluster_dict.has_key(cn):
                        self.cluster_dict[cn]['x'].append(k1)
                        #print k,len(self.line_peak[k]),len(self.cluster_array_dict[k]),num
                        self.cluster_dict[cn]['y'].append(self.line_peak[k1][y1])
                    else:
                        self.cluster_dict[cn] = {'x': [k1],
                                                 'y': [self.line_peak[k1][y1]]
                                                 }

                    x, y = x1, y1
    def peak_cluster(self):
        #best_friend1 = self.best_friend1
        #best_friend2 = self.best_friend2
        line_peak = self.line_peak
        cluster_array = {}
        self.line_number_list = line_peak.keys()
        self.line_number_list.sort()
        for k in self.line_number_list:
            cluster_array[k] = np.zeros(len(line_peak[k]), dtype=int)
        self.cluster_array_dict = cluster_array
        self.cluster_number = 0
        self.peak_cluster_Recursively()

        '''
        ff = open('cluster_dict.pkl','w')
        pickle.dump(self.cluster_dict,ff)
        ff.close()
        sys.exit(0)
        print self.cluster_number
        plt.figure(1,figsize=(20,20))
        ck = self.cluster_dict.keys()
        ck.sort()
        for k in ck:
            x = self.cluster_dict[k]['x']
            y = self.cluster_dict[k]['y']
            plt.plot(y,x,'o-')
        plt.show()
        '''
    def compute_best_friend(self):
        '''
        Generate best_friend1, best_frined2, and Use them to cluster points into fibers
        :return:
        '''
        whole_peaks = self.whole_peaks
        line_peak = self.line_peak
        line_numbers = line_peak.keys()
        line_numbers.sort()
        linen = len(line_numbers)
        best_friend1 = {} # top to bottom
        best_friend2 = {} # bottom to top
        for ind in range(linen-1):
            # i from top to bottom
            # j from bottom to top
            i = line_numbers[ind]
            j = line_numbers[linen - ind -1]
            ii = line_numbers[ind+1]
            jj = line_numbers[linen - ind - 2]

            line1 = line_peak[i]
            line2 = line_peak[ii]
            ret = self.find_best_friend(line1, line2)
            best_friend1[i] = ret

            line3 = line_peak[j]
            line4 = line_peak[jj]
            ret = self.find_best_friend(line3, line4)
            best_friend2[j] = ret
        self.best_friend1 = best_friend1
        self.best_friend2 = best_friend2
        #self.peak_cluster()
        '''
        ind_array = {}

        for ind in range(1,linen-1,1):

            current_line = line_numbers[ind]
            upper_line = line_numbers[ind-1]
            lower_line = line_numbers[ind+1]
            vec = []
            for ii in range(len(line_peak[current_line])):
                upper_best_friend = best_friend2[current_line][ii]
                lower_best_friend = best_friend1[current_line][ii]

                upper_friend_friend = best_friend1[upper_line][upper_best_friend]
                lower_friend_friend = best_friend2[lower_line][lower_best_friend]
                if ii != upper_friend_friend or ii != lower_friend_friend:
                    vec.append(False)
                else:
                    vec.append(True)
            ind_array[current_line] = np.array(vec)
        current_line = line_numbers[0]
        vec = []
        for ii in range(len(line_peak[current_line])):
            lower_line = line_numbers[1]
            lower_best_friend = best_friend1[current_line][ii]
            lower_friend_friend = best_friend2[lower_line][lower_best_friend]
            if ii != lower_friend_friend:
                vec.append(False)
            else:
                vec.append(True)
        ind_array[current_line] = np.array(vec)
        current_line = line_numbers[-1]
        vec = []
        for ii in range(len(line_peak[current_line])):
            upper_line = line_numbers[-2]
            upper_best_friend = best_friend2[current_line][ii]
            upper_friend_friend = best_friend1[upper_line][upper_best_friend]
            if ii != upper_friend_friend:
                vec.append(False)
            else:
                vec.append(True)
        ind_array[current_line] = np.array(vec)

        lp = {}
        wp = []
        for key in line_numbers:

            vec = np.array(line_peak[key])
            vec = vec[ind_array[key]]
            lp[key] = vec
            wp.append((key, lp[key]))
        self.line_peak = lp
        self.whole_peaks = wp
        '''
    def trace_image(self, linelist=None):
        self.find_whole_peaks()
        self.compute_best_friend()
        self.peak_cluster()
        self.find_fibers(polyfit_N=2, threshold_percent=0.4)
        self.fiber = []
        keys = self.detected_fiber.keys()
        keys.sort()
        xx = np.array(range(self.image.image_shape[0]))
        for fib in keys:
            line = self.detected_fiber[fib]
            y_center = np.polyval(line['poly'], xx)
            y = np.array([int(value + 0.5) for value in y_center])
            item = {'x':xx[:], 'y':y, 'y_center':y_center, 'tag':'find'}
            self.fiber.append(item)
        '''
        data = self.image.get_image()
        whole_peaks = self.whole_peaks
        xx = np.array(range(len(data)))

        poly_lines = self.find_trace()
        for pl in poly_lines:
            coeff = pl
            y_center = np.polyval(coeff, xx)
            yy = np.array([int(value+0.5) for value in y_center])
            #plt.plot(yy,xx,'r')
            line = {'x':xx[:],'y':yy,'y_center':y_center, 'tag':'find'}
            self.fiber.append(line)
        '''
        #self.plot_fitting(b, xx, yy, lines)
        #plot_sigma_hist(b, xx, yy, lines)
        #plt.scatter(cx,cy,c='r',edgecolors='none',s=1)
        #plt.show()

    def find_fibers(self, polyfit_N = 2, threshold_percent = 0.4):

        #data = pickle.load(open('./cluster_dict.pkl'))
        data = self.cluster_dict

        #polyfit_N = 2
        #threshold_percent = 0.4

        #print len(data)
        keys = data.keys()
        keys.sort()
        longest_num = 0
        for key in keys:
            if len(data[key]['x']) > longest_num:
                longest_num = len(data[key]['x'])

        #plot_data(data, longest_num)

        detected_fiber = {}


        for key in keys:
            y = data[key]['y']
            if len(y) == longest_num:  # perfect fiber trace
                ind = int(y[0])
                coeff = np.polyfit(data[key]['x'], data[key]['y'],polyfit_N)
                detected_fiber[ind] = {'data':data[key], 'poly':coeff}

                
                data[key]['status'] = 'good'
                data[key]['detected_fiber_key'] = ind
            else:
                data[key]['status'] = 'not good'
                data[key]['detected_fiber_key'] = -1
        #print len(detected_fiber.keys())
        detected_fiber_key_list = detected_fiber.keys()
        detected_fiber_key_list.sort()
        detected_fiber_first_key = detected_fiber_key_list[0]
        line_array = detected_fiber[detected_fiber_first_key]['data']['x']
        not_good_fiber = {}
        not_good_fiber_coeff = {}
        for key in keys:
            if data[key]['status'] == 'good':
                continue
            first_point = (data[key]['x'][0], data[key]['y'][0])
            x_ind = 0
            while x_ind < len(line_array) and first_point[0]>line_array[x_ind]:
                x_ind +=1
            ind = 0
            while ind < len(detected_fiber_key_list) and first_point[1] > detected_fiber[detected_fiber_key_list[ind]]['data']['y'][x_ind]:
                ind +=1
            if ind == len(detected_fiber_key_list):
                kk = int(detected_fiber[detected_fiber_key_list[-1]]['data']['y'][0]+2)
            else:
                kk = int(detected_fiber[detected_fiber_key_list[ind]]['data']['y'][0]-2)
            coeff = []
            if ind >0:
                coeff.append(detected_fiber[detected_fiber_key_list[ind-1]]['poly'])
            else:
                coeff.append(None)
            if ind < len(detected_fiber_key_list):
                coeff.append(detected_fiber[detected_fiber_key_list[ind]]['poly'])
            else:
                coeff.append(None) 
            not_good_fiber_coeff[kk] = coeff
            if not_good_fiber.has_key(kk):
                not_good_fiber[kk].append(key)
            else:
                not_good_fiber[kk] = [key]
        #self.not_good_fiber = not_good_fiber

        #plot_data2(data,detected_fiber, not_good_fiber)

        percent = np.linspace(0.01, 0.99, 100)
        ## possible pairs
        xx_range = detected_fiber[detected_fiber_key_list[0]]['data']['x']
        for key in not_good_fiber.keys():
            rec_dict = {}
            ll = not_good_fiber[key]
            if not_good_fiber_coeff[key][0] is None or not_good_fiber_coeff[key][1] is None:
                continue
            #print key,len(ll)
            points = {'x':[],'y':[]}
            for l in ll:
                points['x'] += list(data[l]['x'])
                points['y'] += list(data[l]['y'])
            points['x'] = np.array(points['x'])
            points['y'] = np.array(points['y'])

            coeff1,coeff2 = not_good_fiber_coeff[key]
            used = []
            ind_list = np.array(range(len(points['x'])))
            number_threshold = int(longest_num * threshold_percent)
            line_list = []
            while True:
                best_percent = -1
                best_num = -1
                best_points_set = []
                for pp in percent:
                    coeff = coeff1 * pp + coeff2*(1-pp)
                    yy = np.polyval(coeff, points['x'])
                    residual = yy - points['y']
                    ind =np.abs(residual) < 1.0
                    near_point_ind = ind_list[ind]
                    near_point_ind = list(near_point_ind)
                    #print near_point_ind
                    i = 0
                    while i<len(near_point_ind):
                        if near_point_ind[i] in used:
                            del near_point_ind[i]
                        else:
                            i+=1
                    num = len(near_point_ind)
                    if num < number_threshold:
                        continue
                    if num < best_num:
                        continue
                    best_num = num
                    best_percent = pp
                    best_points_set = near_point_ind[:]
                if best_num == -1:
                    break   #no good percent
                x,y = [],[]
                for ind in best_points_set:
                    x.append(points['x'][ind])
                    y.append(points['y'][ind])
                x = np.array(x)
                y = np.array(y)
                coeff_new = np.polyfit(x,y,polyfit_N)
                line_list.append((best_num, coeff_new, (x[:],y[:])))
                used += best_points_set
            if len(line_list) == 0:
                continue
            for line in line_list:
                y = np.polyval(line[1], xx_range)
                tmp_d = {'x':line[2][0], 'y':line[2][1]}
                ind = int(y[0])
                detected_fiber[ind] = {'data':tmp_d, 'poly':line[1]}

                #plt.plot(y, xx_range,'g')

        self.detected_fiber = detected_fiber

    '''
    def plot_fitting(self,b, xx, yy, lines, fig=0):
        plt.figure(figsize=(20,5),tight_layout=True)
        plt.clf()
        plt.plot(b)
        plt.plot(xx,yy,'ro')
        for i in range(len(lines)):
            p, gauss_x, gauss_y = lines[i]
            if i%20!=0:
                plt.plot(p[1], np.max(gauss_y),'go')
            else:
                plt.plot(p[1],np.max(gauss_y),'yo')
            plt.plot(gauss_x, gauss_y,'m', alpha=0.5)
        axis = plt.axis()
        vskip = (axis[3]-axis[2])*0.05
        for i in range(len(lines)):
            if i % 20 !=0:
                continue
            p, gauss_x, gauss_y = lines[i]
            plt.text(p[1],np.max(gauss_y)+vskip, '%d'%(i),ha='center')
    '''
    def find_whole_peaks(self, linelist=None):
        '''
        Find all peaks from the image. The peak means the point's value is larger
        than the values of both the left and the right point.
        :param linelist: sampling row number list. if linelist is None, using the sampling_step
        :return:
            self.whole_peaks: all peaks we found. List, cell structral: (line_number, centers)
            self.line_peak: all peaks we found. Dict, {line_number:centers}
        '''
        fiber_half_width = self.fiber_half_width
        cx,cy = [],[]
        whole_peaks = []
        data = self.image.get_image()
        if not isinstance(linelist, list):
            linelist = range(1, len(data),self.sampling_step)
        line_peak = {}
        for i in linelist:
            #i = 981
            lineno = i
            b = data[lineno,:]
            xx,yy = self.find_peak(b)
            lines = self.gaussian_fit(b, xx,yy)
            centers = [ll[0][1] for ll in lines]
            whole_peaks.append((i, centers))
            line_peak[lineno] = centers
        self.whole_peaks = whole_peaks
        self.line_peak = line_peak
    def find_peak(self, y ):
        '''
        Find peaks in a row
        :param y:
        :return:
        '''
        '''
        :param y:
        :return:
        '''
        fiber_half_width=self.fiber_half_width
        x = np.array(range(len(y)))
        gap = fiber_half_width 
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
        return np.array(xx), np.array(yy)
    def gaussian_fit_single(self,x, y, b=None):
        #if b!=801:
        #    return None,None,None,None,-1
        #plt.plot(x,y,'b')
        


        err_ret = (None, None, None,None,-1)

        #if len(y) %2 !=1:
        #    # length of x or y must be an odd number.
        #    return err_ret
        n = len(y) 
        if b == None:
            b = np.mean(x)
        #else:
        #    # The middle one must be the heightest point in the data
        #    tot = np.sum(y<y[int(n/2)])
        #    if tot != n-1:
        #        return err_ret

        gauss_func = lambda x,p: p[0]*np.exp(-0.5*(x-p[1])**2/p[2]**2)+p[3]
        err_func = lambda p, x, y: gauss_func(x,p) - y
        p0 = [np.max(y), b, 3.5,0.0]
        pfit, pcov, infodict, errmsg, success = optimize.leastsq(err_func, p0, args=(x, y), full_output=1)
        
        if (len(y) > len(p0)) and pcov is not None:
            s_sq = (err_func(pfit, x, y)**2).sum()/(len(y)-len(p0))
            pcov = pcov * s_sq
        else:
            pcov = np.inf

        error = [] 
        for i in range(len(pfit)):
            try:
              error.append( np.absolute(pcov[i][i])**0.5)
            except:
              error.append( 0.00 )
        pfit_leastsq = pfit
        perr_leastsq = np.abs(np.array(error))
        
        pfit_leastsq[2] = np.abs(pfit_leastsq[2])

        new_y = gauss_func(x, pfit)
        #plt.plot(x,new_y,'r')
        #print 801,'#############',pfit_leastsq
        return x, new_y, pfit_leastsq, perr_leastsq, success
    def gaussian_fit(self, data, x, y):
        lines = []
        half_width = self.fiber_half_width
        ew_list = []
        for i in range(len(x)):
            
            center_x = x[i]
            center_y = y[i]
            if center_x <half_width or len(data)-center_x<half_width+1:
                #out of the boundary
                continue
            fit_x = np.array(range(center_x-half_width, center_x+half_width+1,1))
            fit_y = data[center_x-half_width:center_x+half_width+1]
            ew = self.compute_ew(fit_x, fit_y)
           
            peak_x = fit_x[half_width]
            g_x, g_y, pfit, perr, success = self.gaussian_fit_single(fit_x, fit_y, peak_x)

            pmag, pcenter, psigma, pbase = pfit

            if success not in [1,2,3,4]:
                # fitting solution not found
                #print_info(center_x, pfit,perr,success)
                continue
            if psigma > 10 or psigma<1:
                # sigma is weired
                #print_info(center_x, pfit,perr,success)
                continue
            #if perr[1] > half_width:
                # fitted center's error is too large
                #print_info(center_x, pfit,perr,success)
            #    continue
            if pcenter <fit_x[0] or pcenter > fit_x[-1]:
                # fitted center is out of range
                #print_info(center_x, pfit,perr,success)
                continue
            if np.abs(pcenter-center_x) > 4:
                # The fitting center is too far away with the peak
                continue
            #if np.abs(np.max(g_y) - np.min(y)) <100.0:
            #if pfit[0] < 300:
                # Signal is too weak
                #print_info(center_x, pfit,perr,success)
            #    continue
            if pfit[0] < 0:
                # Gaussioan function Direction error
                continue
            #if ew<0:
            #    # too low ew is bad
            #    continue
            #print ew
            line = (pfit, g_x, g_y)
            ew_list.append(ew)
            lines.append(line)
        ew_threathold = np.max(ew_list) / 100.0
        i =0
        while i<len(lines):
            ew = ew_list[i]
            line = lines[i]
            if ew <ew_threathold:
                del lines[i]
                del ew_list[i]
                continue
            i+=1
        self.ew_list = ew_list
        return lines
    def compute_ew(self,x,y):
        half_n = int(len(x)/2)
        left_i = y[:half_n].argmin()
        right_i = y[half_n:].argmin() + half_n
        x0, y0 = x[left_i], y[left_i]
        x1, y1 = x[right_i],y[right_i]

            
        #base_line = lambda xx:y0 + 1.0*(xx-x0)*((y1-y0)*1.0/(x1-x0))
        coeff = np.polyfit([x0,x1],[y0,y1],1)
        base_line = np.polyval(coeff, x)
        new_y = y-base_line#(x)

        s = (new_y[1:]+new_y[:-1])*(x[1:]-x[:-1])*0.5
        s = np.sum(s)
        return s
    '''
    def find_trace(self):
        lines = self.whole_peaks
        trace = []
        line_number = 0

        ln = [len(ll[1]) for ll in lines]
        ld = {}
        for number in ln:
            if ld.has_key(number):
                ld[number]+=1
            else:
                ld[number]=1
        line_count_number = 0
        for number in ld.keys():
            if ld[number] > line_count_number:
                line_count_number = ld[number]
                line_number = number
        #print line_number
        new_lines = [ll for ll in lines if len(ll[1]) == line_number]

        self.detected_line_number = line_number
        lines = []
        ##### reshape
        for i in range(line_number):
            x,y = [],[]
            for j in range(len(new_lines)):
                x.append(new_lines[j][0])
                y.append(new_lines[j][1][i])
            lines.append((x,y))
        polylines = []

        for line in lines:
            x,y = line
            coeff = self.robust_polyfit(x,y,5)
            polylines.append(coeff)
        return polylines
    '''
    def robust_polyfit(self,x,y,N):
        coeff = np.polyfit(x,y,N)
        return coeff
        '''
        diff_x, diff_y = [],[]
        i = 0
        while i<len(new_lines)-1:
            for j in range(line_number):
                #diff_x.append(new_lines[i+1][0] - new_lines[i][0])
                diff_x.append(0.5*(new_lines[i+1][0]+new_lines[i][0]))
                diff_y.append(np.abs(new_lines[i+1][1][j]-new_lines[i][1][j]))
            i+=1
        diff_y = np.array(diff_y)
        diff_x = np.array(diff_x)
        ind = diff_y >2.0
        dy = diff_y[ind]
        dx = diff_x[ind]
        plt.scatter(dx,dy,c='r',edgecolors='none',s=3)
        #diff_yy = []
        #i = 0
        #while i<len(diff_y)-1:
        #    diff_yy.append(np.abs(diff_y[i]-diff_y[i+1]))
        #    i+=1
        #plt.hexbin(diff_x,diff_y, mincnt=1)
        
        #plt.colorbar()
        #plt.hist(diff_y)
        #plt.hist(diff_yy)
        '''
        return
    def save(self, ff = None):
        if ff is None:
            return
        if os.path.exists(ff):
            os.remove(ff)
        data = []
        tag = []
        for f in self.fiber:
            data.append(f['y_center'])
            if f['tag']=='find':
                tag.append(0)  #fiber found by trace
            else:
                tag.append(1)  #fiber found by calibration

        hdu = pyfits.PrimaryHDU(np.array(data))
        hdu_tag = pyfits.ImageHDU(np.array(tag))
        hdu_tag.name = 'FIBERTAG'
        hdu_tag.header['TAG_VAL']=('0,1','0:fiber found by trace, 1:by inputed')
        hdulist = pyfits.HDUList([hdu,hdu_tag])
        hdulist.writeto(ff)
    def load(self, ff = None):
        if ff is None:
            return
        if not os.path.exists(ff):
            return
        hdulist = pyfits.open(ff)
        data = hdulist[0].data
        tag = hdulist[1].data
        x = np.array(range(len(data[0])))
        fiber = []
        for i in range(len(data)):
            y_center = data[i]
            y = np.array([int(v+0.5) for v in y_center])
            t = ['find','calibration'][tag[i]]
            item ={'x':x[:], 'y':y, 'y_center':y_center, 'tag':t}
            fiber.append(item)
        self.fiber = fiber


if __name__ == '__main__':
    ff = './tmp/rs-04b-20170414184339-2604-83315203.fit.gz'
    image = image2d.image2d(ff)
    #image.subtract_overscan()
    t = trace(image)
    #t.find_whole_peaks()
    print len(t)

    import matplotlib.pyplot as plt
    plt.figure(0)
    data = image.get_image()
    
    plt.imshow(data, origin='lower')
    '''
    wp = t.whole_peaks
    x,y = [],[]
    for pl in wp:
        for p in pl[1]:
            if isinstance(pl[0],int) and isinstance(p,float):
                x.append(pl[0])
                y.append(p)
            else:
                print pl[0],p
    plt.plot(y,x,'ro')
    plt.figure(1)
    plt.hist(t.ew_list,50)
    '''
    for line in t:
        y = line['x']
        x = line['y']
        plt.plot(x,y,'r')
    plt.show()
    #print t.fiber




