import matplotlib.pyplot as plt
import os
import numpy as np
def plot_trace(t, output='./Plot_trace2d'):
    output_prefix,output_filename = os.path.split(output)
    output_filename = output_filename.replace('.png','')
    image = t.image
    if not image.is_valid():
        return
    image_data = image.get_image()
    
    plt.figure(0,figsize=(50,50))
    plt.clf()
    plt.imshow(image_data, origin='lower')

    wp = t.whole_peaks
    x,y = [],[]
    for row in wp:
        r, centers = row
        for c in centers:
            x.append(r)
            y.append(c)
    plt.plot(y,x,'yo',alpha=0.8)
    n,m = 0,0
    for fib in t:
        x = fib['x']
        y = fib['y_center']
        if fib['tag'] == 'find':
            plt.plot(y,x,'r')
            n+=1
        else:
            plt.plot(y,x,'g')
            m+=1
        if (n+m) %10 == 0:
            plt.text(y[5],x[5],str(n+m),color='w',ha='center',va='center',alpha=.8)
            plt.text(y[-6],x[-6],str(n+m),color='w',ha='center',va='center',alpha=.8)
        
    '''
    detected_fiber = t.detected_fiber
    keys = detected_fiber.keys()
    keys.sort()
    for key in keys:
        lines = detected_fiber[key]['data']
        coeff = detected_fiber[key]['poly']
        # print key, ','.join(['%.3f'%(c) for c in coeff])
        plt.plot(lines['y'], lines['x'], 'ko', alpha=0.2)
        # xx = np.linspace(lines['x'][0],lines['x'][-1],100)
        plt.plot(np.polyval(coeff, lines['x']), lines['x'], 'r')
    '''
    '''
    for key in not_good_fiber.keys():
        ll = not_good_fiber[key]
        x, y = [], []
        for g in ll:
            xx, yy = data[g]['x'], data[g]['y']
            # plt.plot(yy,xx,'o-')

            x += list(xx)
            y += list(yy)
        plt.plot(y, x, 'o')

    wp = t.whole_peaks
    x,y = [],[]
    lines_found = t.detected_line_number
    for pl in wp:
        for p in pl[1]:
            if isinstance(pl[0], int) and isinstance(p, float):
                x.append(pl[0])
                y.append(p)
    plt.plot(y,x,'ro')
    for line in t:
        y = line['x']
        x = line['y_center']
        line_class = line['tag']
        if line_class == 'calibration':
            plt.plot(x,y,'b')
        else:
            plt.plot(x,y,'r')
    
    '''

    plt.title('Lines Detected: %d, Lines Now:%d'%(n, n+m), fontsize = 30)
    plt.savefig(output)

if __name__ == '__main__':
    import image2d
    import trace2d
    ff = './tmp/rs-04b-20170414184339-2604-83315203.fit.gz'
    image = image2d.image2d(ff)
    t = trace2d.trace(image,500)
    plot_trace(t,'/tmp/trace.png')

