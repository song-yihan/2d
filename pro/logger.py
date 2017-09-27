import logging
import os
import sys
import glob
def init(logfile = None, both = True):
    if logfile == None:
        frame = sys._getframe(1) # upper lever function
        fname = frame.f_code.co_filename
        f = os.path.split(fname)[1].replace('.py','.log')
        path = os.getcwd()
        logpath = os.path.join(path,'..','log')
        logfile = os.path.join(logpath,f)
        logfile = logfile.replace('.log','.%d.log'%(os.getpid()))
        logfile_bak = logfile
        while os.path.exists(logfile):
            tmp = logfile +'.*'
            ffs = glob.glob(tmp)
            n = 1
            for ff in ffs:
                ff = ff.replace(logfile,'')[1:]
                if ff == '' or (not ff.isdigit()):
                    continue
                a = int(ff)
                if a+1>n:
                    n=a+1
            logfile = '%s.%d'%(logfile,n)
        if logfile != logfile_bak:
            os.system('mv %s %s'%(logfile_bak, logfile))
            logfile = logfile_bak
        #print "ff:",logpath, logfile
        #print os.getcwd()
        if not os.path.exists(logpath):
            os.mkdir(logpath)
    '''
    #http://blog.csdn.net/jj_liuxin/article/details/3565309
    #logfile = 'log.txt'
             
    logger = logging.getLogger()
    hdlr = logging.FileHandler(logfile)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.NOTSET)
            
    logger.info('message')
    '''
    '''
    print logfile
    #http://blog.chinaunix.net/uid-26000296-id-4372063.html
    logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',  
                    datefmt='%a, %d %b %Y %H:%M:%S',  
                    filename=logfile,  
                    filemode='w')  
    '''
    """Init for logging
    """
    #http://blog.csdn.net/lv26230418/article/details/46356763
    st_format = '%(asctime)s %(filename)s:%(lineno)d %(levelname)s %(message)s'
    st_print = '[%(filename)s:%(lineno)d] %(message)s'
    logging.basicConfig(
                    level    = logging.DEBUG,
                    format   = st_format,
                    datefmt  = '%m-%d %H:%M',
                    filename = logfile,
                    filemode = 'w');
    if not both:
        return
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler();
    console.setLevel(logging.INFO);
    # set a format which is simpler for console use
    formatter = logging.Formatter(st_print);
    # tell the handler to use this format
    console.setFormatter(formatter);
    logging.getLogger('').addHandler(console);

    logging.debug('PROGRAM STARTING')

if __name__ == '__main__':
    print 'CRITICAL > ERROR > WARNING > INFO > DEBUG > NOTSET'
    logger_init()
    logging.debug('test')
    logging.info('infomation')
