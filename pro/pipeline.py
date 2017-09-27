
import sys
import os
import glob

import numpy as np
from scipy import optimize
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import image2d
import trace2d
import extract_1d
import plot2d
import wavelength_calibration

import logger
import logging

class pipeline():
    def __init__(self):
        self.bias = None
        self.bias_list = None
        self.arc = None
        self.arc_list = None
        self.arc_dir_name = None
        pass
    def clear(self):
        pass
    def image_init(self, image_filename_list):
        #if len(image_filename_list) > 5:
        #    image_filename_list = image_filename_list[:5]
        image_list  = [image2d.image2d(ff) for ff in image_filename_list]
        image_list  = [image for image in image_list if image.is_raw_valid()]
        for image in image_list:
            image.subtract_overscan()
        image_list = [image for image in image_list if image.is_valid()]
        return image_list

    def find_trace_ref(self):
        spid_color = self.spid+self.color
        ffs = glob.glob(os.path.join(self.etc_path,'trace_ref','trace_ref_*_%s.fits'%(spid_color)))
        if len(ffs) == 0:
            return None
        ffs.sort()
        dates_st = [os.path.split(ff)[1].split('_')[2] for ff in ffs]
        dates_diff = [np.abs(int(date)-int(dd)) for dd in dates_st]
        dates_closest_index = np.array(dates_diff).argmin()
        return ffs[dates_closest_index]    

    def load_bias(self):
        data_path = self.data_path
        date = self.date
        spid = self.spid
        color =self.color
        bias_path = os.path.join(data_path,date,'bias')
        bias_files = glob.glob(bias_path+'/rb-%s%s-*.fit*'%(spid,color))
        bias_data_list = self.image_init(bias_files)
        if len(bias_data_list) == 0:
            self.bias = None
            self.bias_list = None
            return False
        self.bias = image2d.image_median(bias_data_list)
        self.bias_list = bias_data_list
        return True
    def load_flat(self):
        date_path, date, spid = self.date_path,self.date, self.spid
        color, bias, output_path = self.color, self.bias, self.output_path

        flat_path = os.path.join(data_path,date,'skyflat')
        flat_files = glob.glob(flat_path+'/rs-%s%s-*.fit*'%(spid,color))
        if len(flat_files) == 0:
            self.flat = None
            self.flat_list = None
            return False
        flat_data_list = self.image_init(flat_files)
        flat_data_list = [image-bias for image in flat_data_list]
            
        self.flat = image2d.image_median(flat_data_list)
        self.flat_list = flat_data_list
        logging.info('\nFLAT FILEs:\n\t'+
                 '\n\t'.join([st for st in flat.filename.split(',') if 'bias' not in st]))
        return True

    def write_fits_list(self, image_list, output_subpath, out_filenames = None):
        trace_data = self.trace_data
        if trace_data is None:
            return False
        wave_cali = self.wave_cali

        path = os.path.join(self.output_path, self.date, output_subpath)
        self.make_dir(path)
        ind = 0
        for image in image_list:
            if out_filenames is not None:
                outname = out_filenames[ind]
            else:
                outname = image.filename.split(',')[0]
                outname = os.path.split(outname)[1]
            spec_path = os.path.join(path, outname)
            spec = extract_1d.extract_1d(image)
            spec.set_trace(trace_data)
            if wave_cali is not None:
                spec.set_wavelength(wave_cali)
            spec.extract_all()
            spec.writeto(spec_path)
            del spec
            ind +=1
        return True
    def release_image_list(self, image_list):
        if image_list is None or (not isinstance(image_list, list)):
            return
        for image in image_list:
            image.release_image()
            del image
        
    def input_param(self, **arg):
        for key in arg.keys():
            command = 'self.%s = arg["%s"]'%(key, key)
            exec(command)
    def set_trace(self, trace_file):
        if not os.path.exists(trace_file):
            return
        trace_data = trace2d.trace()
        trace_data.load(trace_file)
        self.trace_data = trace_data
        logging.info('\tUse Trace File:\n%s'%(trace_file))
        logging.info('Traced Fiber Number:%d'%(len(trace_data)))
    def load_trace(self):
        if self.trace_data  is not None:
            return 

        date, spid, color = self.date, self.spid, self.color
        output_path = self.output_path
        trace_exists_file = os.path.join(output_path,date,'trace',
                'trace_data_%s_%s%s.fits'%(date, spid,color))

        if os.path.exists(trace_exists_file) :
            logging.info('\t Using the existed trace file.')
            logging.info('\n%s'%(trace_exists_file))
            trace_file = trace_exists_file
            trace_data = trace2d.trace()
            trace_data.load(trace_file)
            logging.info('Traced Fiber Number:%d'%(len(trace_data)))
            self.trace_data = trace_data
            return 
    def do_trace(self):
        if self.bias is None:
            logging.info('Before Trace, Bias is needed.')
            return
        if self.flat is None:
            logging.info('Before Trace, Flat is needed.')
            return

        flat = self.flat
        trace_data = trace2d.trace(flat)
        logging.info('Traced Fiber Number:%d'%(len(trace_data)))
        
        trace_ref = self.find_trace_ref()
        if trace_ref is not None:
            logging.info('\tTRACE REF:\n%s'%(trace_ref))
            trace_data.calibration(trace_ref)
            logging.info('After Fiber Calibration, Fiber Number:%d'%(len(trace_data)))
        else:
            logging.info('TRACE REF: None')
        png_filename = os.path.join(output_path,date,'trace','trace_%s.png'%(spid+color))
        plot2d.plot_trace(trace_data, png_filename)
        fits_filename = os.path.join(output_path,date,'trace',
                           'trace_data_%s_%s.fits'%(date, spid+color))
        trace_data.save(fits_filename)
        trace_data.remove_image()
        self.trace_data = trace_data

    def init(self):
        self.bias = None
        self.bias_list = None
        self.arc = None
        self.arc_list = None
        self.flat = None
        self.flat_list = None
        self.trace_data = None
        self.wave_cali = None
    def load_arc(self, arc_name = 'arctest1'):

        data_path, date, spid = self.data_path, self.date, self.spid
        color, output_path = self.color,self.output_path
        bias = self.bias
        if self.arc_dir_name is not None:
            arc_name = self.arc_dir_name
        arc_path = os.path.join(data_path, date,arc_name)
        arc_files = glob.glob(arc_path+'/ra-%s-*.fit*'%(spid+color))
        arc_data_list  = self.image_init(arc_files)
        if len(arc_data_list) == 0:
            logging.info('No Arc Files found.')
            arc = None
            self.arc_list = None
            return False
        else:
            arc_data_list = [image-bias for image in arc_data_list]
            self.arc = image2d.image_median(arc_data_list)
            self.arc_list = arc_data_list
        return True
    def do_wavelength_calibration(self):
        #if self.arc_data is None:
        #    logging.info('No Arc Data. Cannot do Wavelength Calibration.')
        #    return
        if self.trace_data is None:
            logging.info('Without trace data, Cannot wavelength calibration')
            return
        arc = self.arc
        if arc is not None:
            logging.info('\nARC FILEs:\n\t'+
                 '\n\t'.join([st for st in arc.filename.split(',') if 'bias' not in st]))
            arc_data = extract_1d.extract_1d(arc)
            arc_data.set_trace(self.trace_data)
            arc_data.extract_all()
            
            logging.info('Wavelength Calibrating....')
            wave_cali = wavelength_calibration.wavelength_calibration(arc_data)
            if color == 'b':
                logging.info('\tUsing Sc Lamp on Blue Wing')
                wave_cali.calibration_sc_blue()
            elif color == 'r':
                logging.info('\tUsing Sc Lamp on Red Wing') 
            logging.info('Wavelength Calibration finished.')
        else:
            wave_cali = None
            logging.info('Load ARC Failed.')
        self.wave_cali = wave_cali
    def fill_tempfile_wavelength(self):
        pass
    def load_wave_cali(self):
        output_path, date = self.output_path, self.date
        spid, color = self.spid, self.color
        wave_cali_file = os.path.join(output_path, date, 'arc', 'wave_cali_%s%s.fits'%(spid, color))
        pass
    def run(self):
        logging.info('='*20+' START '+'='*20)
        logging.info('='*20+'  END  '+'='*20)
    def run_oneday(self):
        date, planid_list, spid,color = self.date, self.planid_list, self.spid, self.color
        data_path,  output_path, etc_path=self.data_path, self.output_path, self.etc_path


        logging.info('Date: %s Spid: %s'%(date, spid))
        self.init()
        logging.info('Loading Bias...')
        self.load_bias()
        if self.bias is None:
            logging.info('Load BIAS Failed.')
            return
        else:
            bias = self.bias
            self.release_image_list(self.bias_list)
            logging.info('\nBIAS FILEs:\n\t'+
                '\n\t'.join(bias.filename.split(',')))
    
        self.load_trace()
        if self.trace_data is None:
            self.load_flat()
            self.do_trace()
            #logging.info('Finding Trace...')
        if self.trace_data is None:
            logging.info('Cannot Trace the flat.')
            return
        
        if self.flat_list is not None:
            self.write_fits_list(self.flat_list, 'flat')
            self.release_image_list(self.flat_list)
            self.flat_list = None
        if self.flat is not None:
            self.write_fits_list([self.flat],'flat','median_flat_%s%s.fits'%(spid,color))
            self.release_image_list([self.flat])
            self.flat = None
        
        self.load_wave_cali()
        if self.wave_cali is None:
            self.load_arc()
            self.do_wavelength_calibration()
            
        if self.arc_list is not None:
            self.write_fits_list(self.arc_list, 'arc')
            self.release_image_list(self.arc_list)
            self.arc_list = None
        if self.arc is not None:
            self.write_fits_list([self.arc],'arc','median_arc_%s%s.fits'%(spid, color))
            self.release_image_list([self.arc])
            self.arc = None

        if self.wave_cali is not None:
            wave_cali_file = os.path.join(output_path, date, 'arc', 'wave_cali_%s%s.fits'%(spid, color))
            #self.wave_cali.writeto(wave_cali_file)
            self.fill_tempfile_wavelength()
            pass
        if not isinstance(planid_list,list):
            logging.info('Planid_list is not a list')
            return
        
        for planid in planid_list:
            logging.info('-'*10+' Planid: '+planid+' '+'-'*10)
            plan_path = os.path.join(output_path, date, planid)
            if not os.path.exists(plan_path):
                os.mkdir(plan_path)
            self.extract_object(planid)


    def extract_object(self,planid):
        data_path, spid, color = self.data_path, self.spid, self.color
        trace_data, wave_cali, output_path=self.trace_data, self.wave_cali,self.output_path
        bias = self.bias
        obj_path = os.path.join(data_path,planid)
        obj_files = glob.glob(os.path.join(data_path,date,planid,'ro-%s-*.fit.gz'%(spid+color)))
        logging.info('\nObj FILEs:\n\t'+'\n\t'.join(obj_files))
        obj_data_list = self.image_init(obj_files)
        obj_data_list = [image-bias for image in obj_data_list]

        for obj in obj_data_list:
            spec = extract_1d.extract_1d(obj)
            spec.set_trace(trace_data)
            if wave_cali is not None:
                spec.set_wavelength(wave_cali)
            spec.extract_all()
            fits_name = obj.filename.split(',')[0]
            fits_name = os.path.split(fits_name)[1]
            fits_name = fits_name.replace('ro-','spec-')
            fits_name = fits_name.replace('fit.gz','fits')
            #print fits_name
            spec.writeto(os.path.join(output_path,date,planid,fits_name))
    def make_filename_dir(self, st):
        d,f = os.path.split(st)
        if d == '' or d == os.path.sep:
            return
        if os.path.exists(d) and os.path.isdir(d):
            return
        else:
            self.make_filename_dir(d)
    def make_dir(self, st):
        d,f = os.path.split(st)
        if d == '' or d == os.path.sep:
            return
        if os.path.exists(d) and os.path.isdir(d):
            return
        else:
            self.make_filename_dir(d)
        os.mkdir(st)

if __name__ == '__main__':
    logger.init()
    

    dates = ['20170609','20170610']
    dates = ['20170610']
    #planid = ['BD47293601','BD26260601']
    spid = '02'
    color = 'b'
    data_path = '/data2/rawdata'
    output_path = '/2dresult'
    etc_path = '/home/syh/work/Project/2d/etc'
    
    pp = pipeline()
    pp.input_param(date='20170610', spid = spid, color= color)
    pp.input_param(data_path = data_path, output_path = output_path)
    pp.input_param(etc_path = etc_path)


    for date in dates:
        path = os.path.join(data_path, date)
        ffs = glob.glob(path+'/*')
        planids = [os.path.split(ff)[1] for ff in ffs if os.path.isdir(ff)]
        planids = [ff for ff in planids if ff not in ['bias','flat','arctest1','arctest2','nightsky','skyflat','arc'] and 'SCHEME' not in ff]
        #print (date, planids)
        #make_output_dir(date, output_path)
        pp.input_param(planid_list = planids)
        pp.run_oneday()

    logging.info('PROGRAM END')
