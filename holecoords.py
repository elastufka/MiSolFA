#EL 06.06.17
#change coordinates of hole measurements to uniform coordinate system
import matplotlib.pyplot as plt

def parse_text(filename,circle=True,points=True,distance=True,lines=False):
    '''Parse the output files to get coordinates of circles and/or windows or whatever else'''
    params={'Item':[],'x':[],'y':[],'z':[],'R':[],'e':[]}
    numlines=sum(1 for line in open(filename, 'rU'))
    with open(filename,'rU') as f:
        for i in range(1,numlines):
            try:
                line=f.next()
            except StopIteration:
                break
            #print line
            #if line.strip() == '': line=f.next()
            if circle and line.startswith('Circle'):
                params['Item'].append(line[line.find('-')+1:line.find('(')].strip())
                #print f.tell(), type(f.tell()),f.seek(4),type(f.seek(f.tell()))
                line=f.next()
                params['x'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['y'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['z'].append(float(line[line.find('=')+1:].strip()))
                f.next()
                line=f.next()
                params['R'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['e'].append(float(line[line.find('=')+1:].strip()))
            if points and line.startswith('Point'):
                params['Item'].append(line[line.find(':')+1:line.find('(')].strip())
                #print f.tell(), type(f.tell()),f.seek(4),type(f.seek(f.tell()))
                line=f.next()
                params['x'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['y'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['z'].append(float(line[line.find('=')+1:].strip()))
                params['R'].append('')
                params['e'].append('')                
            if distance and line.startswith('Distance'):
                params['Item'].append(line[line.find(':')+1:].strip())
                #print f.tell(), type(f.tell()),f.seek(4),type(f.seek(f.tell()))
                line=f.next()
                params['x'].append(float(line[line.find('=')+1:].strip())) #note this is now Dx
                line=f.next()
                params['y'].append(float(line[line.find('=')+1:].strip())) #note this is now Dy
                params['z'].append('')                
                params['R'].append('')
                params['e'].append('')                
            if lines and line.startswith('Line'):
                params['Item'].append(line[line.find('-')+1:line.find('(')].strip())
                #print f.tell(), type(f.tell()),f.seek(4),type(f.seek(f.tell()))
                line=f.next()
                params['x'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['y'].append(float(line[line.find('=')+1:].strip()))
                line=f.next()
                params['z'].append(float(line[line.find('=')+1:].strip()))
                params['R'].append('')
                params['e'].append('')                
                
    print params
    return params

def change_coords(params):
    '''Change coordinate system to one relative to the center of the segments? to make it easier on the machining people?'''

#    foo

def write_csv(params,csvname):
    import csv
    with open(csvname,'wb') as f:
        w=csv.writer(f)
        #w.writerow(params.keys())
        w.writerow(['Item','x','y','z','R','e'])
        for row in zip(params['Item'],params['x'],params['y'],params['z'],params['R'],params['e']):
            w.writerow(row)

    
def plot_all(params,labels=True):
    '''Plot everything so I can determine by eye what's what'''
    fig=plt.figure()
    ax=fig.add_subplot(111)
    #for pair in zip(params['x'],params['y']):

    #filter distances out
    #params={i:params[i] for i in params if params[].startswith('Distance')} #fix this..
    #print params['Item']
    #ax.scat()
    ax.scatter(params['x'],params['y'])
    if labels:
        for i,x,y in zip(params['Item'],params['x'],params['y']):
            ax.annotate('%s' % i,xy=(x,y),textcoords='data')
        plt.grid()
    fig.show()
#if __name__ != main:
#filename='/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02/mw649sub2765_2017_06_02.txt'
#params=parse_text('/Users/wheatley/Documents/Solar/MiSolFA/prototypes/mw469sub2765_2017_06_02/'+filename)
#plot_all(params)
