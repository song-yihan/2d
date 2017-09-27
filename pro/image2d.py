import astropy.io.fits as pyfits
import numpy as np
#import filename_parser
import os
import logging
import jd2gd
import glob
class image2d():
    def __init__(self, filename='', raw_image_shape=(4136,4160), image_shape=(4136,4096),etc_path = '/home/syh/work/Project/2d/etc'):
        if not isinstance(filename, str):
            filename = ''
        self.filename = filename
        if filename != '':
            logging.debug('Loading 2d File:%s'%(filename))
        #self.fp = filename_parser.filename_parser()
        self.raw_image_shape = raw_image_shape
        self.image_shape = image_shape
        self.etc_path = etc_path

        self.raw_image = None #np.zeros(self.raw_image_shape)
        self.image = None #np.zeros(self.image_shape)
        self.error = None #np.zeros(self.image_shape)
        self.invvar = None
        self.header = None

        ########  init parameters, which change when loading fits ##########
        self.obsdate = '20111024' 
        self.left_readoutnoise = 3.0 
        self.right_readoutnoise = 3.0 
        self.left_gain = 1.0  
        self.right_gain = 1.0 
        self.spid = '01'
        self.color = 'b'
        self.mjm = 80436557
        self.mjm_obsdate = '20111024'
        self.filename_obsdate = '20111024'
        self.header_obsdate = '20111024'
        self.mjm_obsdate = '20111024'
        ############## end init ###################

        self.load_raw_image()
    def release_image(self):
        del self.image
        del self.raw_image
        del self.error
        del self.invvar
        del self.header
    def parse_filename(self):
        filename = self.filename

        raw_filename = os.path.split(filename)[1]
        #'rs-05b-20140111164950-2771-81602938.fit.gz'
        items = raw_filename.split('-')
        spid = items[1][:2]
        color = items[1][2]
        filename_obsdate = items[2][:8]
        mjm = int(items[4].split('.')[0])
        mjm_obsdate = jd2gd.jd2gd(mjm/60.0/24.0)[:8]
        
        self.spid = spid
        self.color = color
        self.mjm = mjm
        self.filename_obsdate = filename_obsdate
        self.mjm_obsdate = mjm_obsdate
    def obtain_amp_direction(self):
        self.amp_direction = 'left_right'
    def load_raw_image(self):

        filename = self.filename
        if filename == '':
            return
        if not isinstance(filename,str):
            return
        if not os.path.exists(filename):
            return
        
        try:
            hdulist = pyfits.open(filename)
            image = hdulist[0].data
            header = hdulist[0].header
            
        except:
            logging.debug('Fits file load Error, filename:%s'%(filename))
            return
        
        self.raw_image = image
        self.header = header
        self.image = image[0:self.image_shape[0],0:self.image_shape[1]]
        self.error = self.image * 0.0
        self.invvar = self.image * 0.0

        self.parse_filename()
        self.obtain_amp_direction()

        head = dict(header)
        if head.has_key('DATE-OBS'):
            self.header_obsdate = head['DATE-OBS'][:10]
            self.header_obsdate = self.header_obsdate.replace('-','')
        ##############################################
        # The Gains in the fits head are not used in 
        # LAMOST pipeline 
        #if head.has_key('GAIN1'):
        #    self.fits_gain1 = float(head['GAIN1'])
        #if head.has_key('GAIN2'):
        #    self.fits_gain2 = float(head['GAIN2'])
        #############################################
        if self.header_obsdate != self.mjm_obsdate or self.header_obsdate != self.filename_obsdate:
            #logging.debug('Obsdate in Fits file are not same')
            #return 
            self.obsdate = self.mjm_obsdate
        else:
            self.obsdate = self.mjm_obsdate
        mean_gain = 1.0
        if os.path.exists(self.etc_path):
            # load gain and read out noise
            ccd_info_files = glob.glob(os.path.join(self.etc_path,'ccdinfo_*.txt'))
            ccd_info_dates = [int(ff.split('_')[-1].replace('.txt','')) for ff in ccd_info_files]
            ccd_info_files.sort()
            ccd_info_dates = np.array(ccd_info_dates)
            ind_max_less_array = ccd_info_dates < int(self.obsdate)
            ind_max_less = ind_max_less_array.argsort()[-1]
            ccd_info_file = ccd_info_files[ind_max_less]
            for line in open(ccd_info_file):
                line = line.strip()
                line = line.replace('\t',' ')
                while '  ' in line:
                    line = line.replace('  ',' ')
                items = line.split(' ')
                if len(items)!=5:
                    continue
                if items[0] != self.spid+self.color.upper():
                    continue
                values = [float(item) for item in items[1:]]
                self.left_gain,self.right_gain,self.left_readoutnoise,self.right_readoutnoise = values
                mean_gain = 0.5*(self.left_gain + self.right_gain)
                break

        self.invvar = 1.0/(np.abs(self.image*mean_gain))
        self.error = np.sqrt(np.abs(self.image*mean_gain))
        return


    def get_image(self):
        
        return self.image
    def get_raw_image(self):
        return self.raw_image
    '''
    def load(self, filename = None):

        self.load_fits(filename)
        if self.is_error():
            return
        
        self.check_image()
        if self.is_error():
            return
        
        self.remove_overscan()
        if self.is_error():
            return



        pass
    '''
    def subtract_overscan(self):
        if not self.is_raw_valid():
            logging.debug('Can not subtract overscan, because raw data is not valid. ')
            return
        logging.debug('Entering subtract_overscan')
        logging.debug('Doing Subtract Overscan, amp_left_right')
        self.amp_left_right()
        logging.debug('Leaving subtract_overscan')
    def measure_readnoise(self, vector):
        s = vector.shape
        n = 1
        for v in s:
            n = n*v
        a = vector.reshape(n)
        a.sort()
        percent = int(n*0.1)
        mean_value = np.mean(a[percent:-percent])
        median_value = np.median(a[percent:-percent])
        if np.abs(mean_value-median_value) >= np.abs(median_value)*0.05:
            logging.warning('Read Noise in Overscan seems Wierd.')
            mean_value = median_value
        s1 = np.std(a[percent:-percent])
        index = (np.abs(a)-mean_value)<7.5*s1
        valid_data = a[index]
        noise = np.std(valid_data)
        return noise
    def amp_left_right(self):
        logging.debug('Using amp_left_right')
        '''
            4135  XXXXXXXXYYYYYYYYAAAAABBBBB
                  XXXXXXXXYYYYYYYYAAAAABBBBB
                  XXXXXXXXYYYYYYYYAAAAABBBBB
                  XXXXXXXXYYYYYYYYAAAAABBBBB
                  XXXXXXXXYYYYYYYYAAAAABBBBB
                  XXXXXXXXYYYYYYYYAAAAABBBBB
                  XXXXXXXXYYYYYYYYAAAAABBBBB
              0   XXXXXXXXYYYYYYYYAAAAABBBBB
                  0     2048    4096 4128 4159
            
            XX is the Left part Data. YY is the Right part Data.
            The length of A and B are the same.

        '''
        #print self.data['image'].shape
        image = self.raw_image
        data = np.zeros(self.image_shape)
        self.image = data
        self.left_readoutnoise = self.measure_readnoise(image[100:4035,4096:4128])
        self.right_readoutnoise = self.measure_readnoise(image[100:4035,4128:4159])
        #shape = self.raw_image_shape
        length_overcan = int((4159-4095)/2.0)
        data_left      = np.array([True]*2048+[False]*2112)
        data_right     = np.array([False]*2048+[True]*2048+[False]*64)
        overscan_left  = np.array([False]*4096+[True]*32+[False]*32)
        overscan_right = np.array([False]*4096+[False]*32+[True]*32) 
        #s = self.data['image'].shape
        for i in range(self.image_shape[0]):
            '''
                For Left Part:
            '''
            dleft = image[i][data_left] * self.left_gain
            oleft = image[i][overscan_left]
            v = self.vector_subtract_overscan(dleft, oleft)
            data[i][data_left[:4096]] = v
            errvalue = self.left_readoutnoise + np.abs(v)
            self.error[i][data_left[:4096]] = np.sqrt(errvalue)
            self.invvar[i][data_left[:4096]] = 1.0/errvalue
            '''
                For Right Part:
            '''
            dright = image[i][data_right] * self.right_gain
            oright = image[i][overscan_right]
            v = self.vector_subtract_overscan(dright, oright)
            data[i][data_right[:4096]] = v
            errvalue = self.right_readoutnoise + np.abs(v)
            self.error[i][data_right[:4096]] = np.sqrt(errvalue)
            self.invvar[i][data_right[:4096]] = 1.0/errvalue
        #print np.sum(data)
        #self.data['data'] = data
        pass

    def amp_top_bottom(self):
        pass
    def vector_subtract_overscan(self, vector, overscan):
        '''
            vector and overscan are both vector.
            call self.compute_overscan(overscan) to compute a reasonable value.
            vector minus the value to get the data
        '''
        overscan_value, err = self.compute_overscan_value(overscan)

        vector = vector - overscan_value
        return vector
    def compute_overscan_value(self, overscan):
        n = len(overscan)
        mean = np.median(overscan) # median replace the mean
        std = np.std(overscan)
        err = 0
        return mean, err
    def is_raw_valid(self):
        if self.raw_image is None:
            return False
        if not isinstance(self.raw_image, np.ndarray):
            return False
        if self.raw_image.shape != self.raw_image_shape:
            return False
        return True
        
    def is_valid(self):
        if self.image is None:
            return False
        if not isinstance(self.image, np.ndarray):
            return False
        if self.image.shape != self.image_shape:
            return False
        return True
    #def __radd__(self, other):
    #    return None
    def __add__(self,other):
        if isinstance(other, image2d):
            if not self.is_valid() or not other.is_valid():
                return None
            ret = image2d()
            ret.filename = ','.join(self.filename, other.filename)
            ret.image = self.image + other.image
            error_value = self.error**2 + other.error**2
            ret.error = np.sqrt(error_value)
            ret.invvar = 1.0/ error_value
            return ret
        return self
    def __sub__(self,other):
        if isinstance(other, image2d):
            if not self.is_valid() or not other.is_valid():
                return None
            ret = image2d()
            ret.filename = ','.join([self.filename, other.filename])
            ret.image = self.image - other.image
            error_value = self.error**2 + other.error**2
            ret.error = np.sqrt(error_value)
            ret.invvar = 1.0/ error_value

            return ret
        return self
    def __div__(self, other):
        if isinstance(other, int):
            other = other * 1.0
        if isinstance(other, float):
            if not self.is_valid() :
                return None
            ret = image2d()
            ret.filename = self.filename
            ret.image = self.image / other

            return ret
        return self
    def __str__(self):
        l = []
        l.append('')
        l.append('2d file name is %s'%(self.filename))
        l.append('2d raw data is valid: %s'%(str(self.is_raw_valid())))
        l.append('2d data(after sub-overscan) is valid: %s'%(str(self.is_valid())))
        l.append('The shape of raw data is: %s'%(str(self.raw_image_shape)))
        l.append('The shape of data is: %s'%(str(self.image_shape)))
        l.append('')
        st = '\n'.join(l) 
        return st
    def clean(self):
        pass
def image_median(image_list = None):
    if not isinstance(image_list, list):
        return image2d()
    image_list = [image for image in image_list if isinstance(image, image2d) and image.is_valid()]
    if len(image_list)==0:
        return image2d()

    image = image2d()

    image.filename = ','.join([img.filename for img in image_list])
    #for k in ['error','info','warning']:
    #    ret.message[k] = self.message[k] + other.message[k]
    
    '''
    shape = image_list[0].image_shape
    image.image = np.zeros(shape)
    for i in range(shape[0]):
        for j in range(shape[1]):
            vector = [img.image[i,j] for img in image_list]
            value = np.median(vector)
            image.image[i,j] = value
    '''
    tmp = np.array([item.image for item in image_list])
    image.image = np.median(tmp, axis=0)
    error = image_list[0].error * 0
    for aa in image_list:
        error = error + aa.error**2
    image.error = np.sqrt(error)
    image.invvar = 1.0 / error
    return image

        
if __name__ == '__main__':
    ff = '../data/20170414/data2/skyflat/rs-04b-20170414184339-2604-83315203.fit.gz'
    image = image2d(ff)
    
    print image
    image.subtract_overscan()
    print image

    '''
    import matplotlib.pyplot as plt

    plt.figure(0)
    plt.imshow(d, origin='lower')
    plt.figure(1)
    print d.shape
    a = d.reshape(4096*4136)
    plt.hist(a,50)
    plt.show()
    '''
