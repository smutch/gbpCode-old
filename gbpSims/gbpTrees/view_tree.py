#!/usr/bin/env python2.6
from numpy import *
import matplotlib
#matplotlib.use('Pdf')
#matplotlib.use('ps')
import matplotlib.pyplot as P
from matplotlib.pyplot import figure, show
from matplotlib.image import NonUniformImage
from matplotlib import colors, ticker
from matplotlib.mlab import bivariate_normal
import pylab 

######################
class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
 
  Register this function like this:
   
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, pos=5, axis=None, xtol=None, ytol=None):
   self.data = zip(xdata, ydata, annotes)
   if xtol is None: 
     xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
   if ytol is None: 
     ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
   self.xtol = xtol
   self.ytol = ytol
   self.pos  = pos
   if axis is None: 
     self.axis = pylab.gca()
   else:
     self.axis= axis
   self.drawnAnnotations = {}
   self.links = []

  def distance(self, x1, x2, y1, y2):
    """
    return the distance between two points
    """
    return math.hypot(x1 - x2, y1 - y2)

  def __call__(self, event):
   if event.inaxes:
     clickX = event.xdata
     clickY = event.ydata
     if self.axis is None or self.axis==event.inaxes:
       annotes = []
       for x,y,a in self.data:
         if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
           annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
       if annotes:
         annotes.sort()
         distance, x, y, annote = annotes[0]
         self.drawAnnote(event.inaxes, x, y, annote)
         for l in self.links:
           l.drawSpecificAnnote(annote)

  def drawAnnote(self, axis, x, y, annote):
   """
   Draw the annotation on the plot
   """
   if (x,y) in self.drawnAnnotations:
     markers = self.drawnAnnotations[(x,y)]
     for m in markers:
       m.set_visible(not m.get_visible())
     self.axis.figure.canvas.draw()
   else:
     print annote
     if(self.pos==1):
       t = axis.text(x,y,annote+' ',ha='right',va='bottom',fontsize=10)
     elif(self.pos==2):
       t = axis.text(x,y,annote,ha='center',va='bottom',fontsize=10)
     elif(self.pos==3):
       t = axis.text(x,y,' '+annote,ha='left',va='bottom',fontsize=10)
     elif(self.pos==4):
       t = axis.text(x,y,annote+' ',ha='right',va='center',fontsize=10)
     elif(self.pos==5):
       t = axis.text(x,y,annote,ha='center',va='center',fontsize=10)
     elif(self.pos==6):
       t = axis.text(x,y,' '+annote,ha='left',va='center',fontsize=10)
     elif(self.pos==7):
       t = axis.text(x,y,annote+' ',ha='right',va='bottom',fontsize=10)
     elif(self.pos==8):
       t = axis.text(x,y,annote,ha='center',va='bottom',fontsize=10)
     elif(self.pos==9):
       t = axis.text(x,y,' '+annote,ha='left',va='bottom',fontsize=10)
     else:
       t = axis.text(x+self.xtol,y,annote,ha='center',va='center')
     m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
     self.drawnAnnotations[(x,y)] =(t,m)

     self.axis.set_autoscale_on(False) 
     self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
   annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
   for x,y,a in annotesToDraw:
     self.drawAnnote(self.axis, x, y, a)
########### End AnnoteFinder class ###########

# Define the tree structure
tree = [('descendant',      'i'),
        ('progenitor_first','i'),
        ('progenitor_next', 'i'),
        ('group_halo_first','i'),
        ('group_halo_next', 'i'),
        ('n_particles',     'i'),
        ('mass',            'f'),
        ('r_vir',           'f'),
        ('x',               'f'),
        ('y',               'f'),
        ('z',               'f'),
        ('vx',              'f'),
        ('vy',              'f'),
        ('vz',              'f'),
        ('sigma_v',         'f'),
        ('v_max',           'f'),
        ('spin_x',          'f'),
        ('spin_y',          'f'),
        ('spin_z',          'f'),
        ('id_MBP',          'int64'),
        ('junk',            'i'),
        ('snap_num',        'i'),
        ('halo_index',      'i'),
        ('halo_id',         'i'),
        ('group_id',        'i')]

tree_mil = [('descendant',      'i'),
            ('progenitor_first','i'),
            ('progenitor_next', 'i'),
            ('group_halo_first','i'),
            ('group_halo_next', 'i'),
            ('n_particles',     'i'),
            ('mass',            'f'),
            ('mass2',           'f'),
            ('mass3',           'f'),
            ('x',               'f'),
            ('y',               'f'),
            ('z',               'f'),
            ('vx',              'f'),
            ('vy',              'f'),
            ('vz',              'f'),
            ('sigma_v',         'f'),
            ('v_max',           'f'),
            ('spin_x',          'f'),
            ('spin_y',          'f'),
            ('spin_z',          'f'),
            ('id_MBP',          int64),
            ('snap_num',        'i'),
            ('file_num',        'i'),
            ('halo_index',      'i'),
            ('half_mass',       'f')]

def tree_spacing():
  return 20
###

def point_size(n_particles):
  return 2.5*log10(n_particles)+2
###

def point(x,y,annotation,r,p,c,cmap,ytol):
  a=[x]
  b=[y]
  n=[annotation]
  s=[r*r]
  color=matplotlib.colors.colorConverter.to_rgb(cmap(1.+log10(c))) 
  #sctpt=ax.scatter(a,b,s=s,color=color,cmap=cmap,edgecolors='black')
  sctpt=pylab.gca().scatter(a,b,s=s,color=color,cmap=cmap,alpha=1.,picker=5)
  af1 = AnnoteFinder(a,b,n,xtol=r,ytol=ytol,pos=p)
  P.connect('button_press_event',af1)
  return sctpt
###

def branch_range(halos,index,offset,i_prog):
  n_particles=halos[index]['n_particles']
  progenitor =halos[index]['progenitor_first']
  offset_i   =offset
  j_prog     =0
  while progenitor>=0:
    t_range   =branch_range(halos,progenitor,offset_i,j_prog)
    offset_i  =t_range[1]
    j_prog    =j_prog+1
    progenitor=halos[progenitor]['progenitor_next']
  if(j_prog==0):
    offset_i=offset_i+tree_spacing()
  return offset,offset_i
###

def plot_tree(halos,index,offset,a_list,n_particles_desc,cmap):
  descendant =halos[index]['descendant']
  halo_index =halos[index]['halo_index']
  snap_num   =halos[index]['snap_num']
  n_particles=halos[index]['n_particles']
  n_snap     =len(a_list)
  delta_a    =1./n_snap

  # Determine plot-position of descendants
  #   (needed before plotting progenitors so we can plot lines)
  t_range=branch_range(halos,index,offset,0)
  x_desc =(t_range[0]+t_range[1])/2.
  y_desc =a_list[snap_num]

  # Plot lines
  progenitor=halos[index]['progenitor_first']
  offset_i  =offset
  snap_min  =snap_num
  snap_max  =snap_num
  count     =0
  halo_count=0
  while progenitor>=0:
    # Walk tree
    pos_prog=plot_tree(halos,progenitor,offset_i,a_list,n_particles,cmap)

    # Set new offset
    offset_i=pos_prog[6]

    # Keep track of max/min snapshot numbers
    snap_min=min([snap_min,pos_prog[3]])
    snap_max=max([snap_max,pos_prog[4]])

    # Count halos under this node
    halo_count=halo_count+pos_prog[5]

    # Plot lines
    x=[pos_prog[0],x_desc]
    y=[pos_prog[1],y_desc]
    P.plot(x,y,'black')

    count=count+1
    progenitor=halos[progenitor]['progenitor_next']

  # Plot points
  colour=float(n_particles)/n_particles_desc
  if(colour>1.):
    colour=1.
  if(colour<0):
    colour=0.
  annotation=' ('+str(index)+'->'+str(descendant)+')('+str(halo_index)+'/'+str(snap_num)+')\n  $n_p$='+str(n_particles)+'('+str(halos[index]['group_halo_first'])+'/'+str(halos[index]['group_halo_next'])+')'
  annotation=' ('+str(a_list[snap_num])+'/'+str(halos[index]['halo_id'])+'/'+str(a_list[halos[descendant]['snap_num']])+str(halos[descendant]['halo_id'])+')('+str(halos[index]['group_id'])+'->'+str(halos[descendant]['group_id'])+')('+str(halos[index]['mass'])+'/'+str(halos[index]['r_vir'])+'/'+str(halos[index]['v_max'])+'/'+str(halos[index]['x'])+'/'+str(halos[index]['y'])+'/'+str(halos[index]['z'])+'/'+str(halos[index]['vx'])+'/'+str(halos[index]['vy'])+'/'+str(halos[index]['vz'])+')'

  point(x_desc,y_desc,annotation,point_size(n_particles),6,colour,cmap,0.5*delta_a)

  # Add this node to the halo count
  halo_count=halo_count+1

  return x_desc,y_desc,x_desc,a_list[snap_min],a_list[snap_max],halo_count,t_range[1]
###

########### MAIN ###########
#filename   = 'trees_930.0'
#filename_a = 'trees_930.0.a_list'
filename   = '/nfs/dset/shrek071/millenium/bolshoi_joinon/treedata/trees_175.2'
filename_a = '/nfs/dset/shrek071/millenium/bolshoi/treedata/Bolshoi.a_list'
tree_dtype = dtype(tree)
#filename   = '/nfs/dset/shrek040/millenium/simulations/millennium_mini/treedata/trees_063.0'
#filename_a = '/nfs/cluster/evo/gpoole/src/3rd_Party/SAM/input/millennium_zlist.txt'
#tree_dtype = dtype(tree_mil)

tree_read  = 4
cmap       = P.matplotlib.cm.jet

tree_read=tree_read+1

# Read list of snapshot expansion factors
a_list=loadtxt(filename_a,'f')

# Open file and read header information
fd          =open(filename,'rb')
n_trees     =fromfile(file=fd,dtype='i',count=1)
n_halos     =fromfile(file=fd,dtype='i',count=1)
n_halos_tree=fromfile(file=fd,dtype='i',count=n_trees)

print 'file            =',filename
print 'selected_tree   =',tree_read
print 'n_halos         =',n_halos
print 'n_trees         =',n_trees
print
print 'n_halos_selected=',n_halos_tree[tree_read-1]

# Skip unwanted trees
for i in range(0,tree_read-1):
	for j in range(0,n_halos_tree[i]):
		halo=fromfile(file=fd,dtype=tree_dtype,count=1)

# Read the desired tree
halos=fromfile(file=fd,dtype=tree_dtype,count=n_halos_tree[tree_read-1])

# Count the number of subtrees
n_subtrees=0
for i in range(0,n_halos_tree[tree_read-1]):
#  print halos[i]
  descendant=halos[i]['descendant']
  print ' ('+str(a_list[halos[i]['snap_num']])+'/'+str(halos[i]['halo_id'])+'/'+str(a_list[halos[descendant]['snap_num']])+'/'+str(halos[descendant]['halo_id'])+')('+str(halos[i]['group_id'])+'->'+str(halos[descendant]['group_id'])+')('+str(halos[i]['mass'])+'/'+str(halos[i]['r_vir'])+'/'+str(halos[i]['v_max'])+'/'+str(halos[i]['x'])+'/'+str(halos[i]['y'])+'/'+str(halos[i]['z'])+'/'+str(halos[i]['vx'])+'/'+str(halos[i]['vy'])+'/'+str(halos[i]['vz'])+')'
  if(halos[i]['descendant']==-1):
    n_subtrees=n_subtrees+1
print 'n_subtrees      =',n_subtrees
# Initialize the plot
fig=figure()
ax=fig.add_subplot(111,axisbg='white')
ax.set_position([0.1,0.07,0.72,0.86])
ax.xaxis.set_visible(False)
ax.set_autoscale_on(False)
ax.set_ylabel('Expansion Factor')

# Plot trees
halo_start=0
halo_count=0
offset    =0.
i_subtree =0
while(halo_start<n_halos_tree[tree_read-1]):
  while(halos[halo_start]['descendant']!=-1 and halo_start<n_halos_tree[tree_read-1]-1):
    halo_start=halo_start+1
  if(halos[halo_start]['descendant']==-1):
    rval=plot_tree(halos,halo_start,offset,a_list,halos[halo_start]['n_particles'],cmap)
    halo_count+=rval[5]
    offset     =rval[6]
    i_subtree  =i_subtree+1
    print rval[5],'halos from subtree',i_subtree,'have been plotted (',halo_count,'in total so far; start=',halo_start,')'
  halo_start=halo_start+1

# Report number of halos plotted
print 'A total of',halo_start,'halos plotted'
print 'Done.'

# Plot a colourbar
a =linspace(1,0.1,256).reshape(-1,1)
ay=fig.add_subplot(1,13,13)
im=ay.imshow(a,aspect='auto',cmap=cmap,extent=[0,1,0.1,1.0])
ay.axes.get_xaxis().set_visible(False)
ay.axes.get_yaxis().tick_right()
ay.axes.get_yaxis().set_label_position("right") 
ay.set_ylabel('$log_{10}$(Merger Ratio [$m_{prog}/m_{desc}$])',fontsize=15)

# Set initial plot limits
ax.set_xlim(-0.05*offset,1.05*offset)
ax.set_ylim(0.,1.05)

ax.set_autoscale_on(False)
ay.set_autoscale_on(False)

fig.savefig('tree.pdf')
#show()
