from PIL import Image

def get_coords(filename):
    x,y,z=[],[],[]
    #pixsize=2.5 #is this right? objective is 2.5x
    with open('mosaic_coords.txt') as f:
    #with open(filename) as f:
        lines=f.readlines()
        for l in lines:
            x.append(float(l[2:9]))
            y.append(float(l[12:18]))
            z.append(float(l[21:27])) #probably need to center these too relative to one coord
        return x,y,z

def fix_coords(x,y,z,ref): #ref is coordinates of top-left corner [198.93.5]
    pixsize=2.5 #is this right? objective is 2.5x -> 1mm -> 2.5 mm (don't think this is right)
    #image size is 640x480
    x0=ref[0]
    y0=ref[1]
    for i in range(0,len(x)):
        print 'before',x[i],y[i]
        x[i]=((x[i]-x0)/pixsize)*640.
        y[i]=((y0-y[i])/pixsize)*480.
        print 'after',x[i],y[i]
    #x[ref]=0.
    #y[ref]=0.
    return x,y,z
        
#https://bytes.com/topic/python/answers/22566-overlaying-transparent-images-pil        
def make_mosaic(x,y,z,filenames=None):
    if not filenames:
        import glob
        filenames=glob.glob('mosaic_h03/*.tif')

    background=Image.new("1",[1408,1152], "white") #some image of some size
    #box1=[640,480,1280,960]
    #background.paste(firstim,box1)
    for i,fn in enumerate(filenames):
        #read in the image - convert them to png first using convert
        im1=Image.open(fn)
        #box=(im1.size[0]-int(x[i]),im1.size[1]-int(y[i]),2*im1.size[0]-int(x[i]),2*im1.size[1]-int(y[i]))
        box=(int(x[i]),int(y[i]),im1.size[0]+int(x[i]),im1.size[1]+int(y[i]))
        print fn, box
        background.paste(im1,box)
    i=3
    im1=Image.open(filenames[i])
    box=(int(x[i]),int(y[i]),im1.size[0]+int(x[i]),im1.size[1]+int(y[i]))
    background.paste(im1,box)
    
    background.show()
    background.save('h03-mosaic.bmp')

if __name__ != "main":
    a,b,c=get_coords('mosaic_coords.txt')
    xx,yy,zz=fix_coords(a,b,c,[198.,93.5])
    foo=make_mosaic(xx,yy,zz)
