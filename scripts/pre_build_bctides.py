#!/home/khan/.Anaconda2/bin/python
# -*- coding=utf-8
#
#       Compute bctides.in à partir des données
#	    FES2012	 
#
#       Gael Arnaud - 12 janvier 2017 
#  
# 2017/07/12 - modifications Yann Krien
# 2017/07/18 - cleaning up Jamal Khan #
# 2023/12/05 - update to python 3 Jamal Khan
########################################################################

from pylab import *
from netCDF4 import Dataset
from pyproj import Proj, transform
import os
import numpy as np

hgrid_file_name = 'hgrid.gr3'
bctide_example = 'bctides.in.sample'
bctide_outfile = 'bctides.ini'
global_atlas_dir = '/run/media/khan/Workbench/Data/FES2012/'

l = os.listdir(global_atlas_dir)
path = global_atlas_dir

class BcTides(object):
  """
  Class bctides
	bctides.in file with 
	48 characters date
	List of Tidal Potential Constituants
	List of Boundary Forcing Frenquencies
	Number of Open Boudaries
  """
  def __init__(self,date=None,tip_dp=None,list_tip=None,list_bfr=None,list_openBound=None):
    self.date=date
    self.tip_dp=tip_dp
    self.list_tip=list_tip
    self.list_bfr=list_bfr 
    self.list_openBound=list_openBound
    #self.flags=flags

class BcFlag(object):
  """
  Class of flags for BcTides class 
  """
  def __init__(self,neta=None,iettype=None,ifltype=None,itetype=None,isatype=None,itrtype=' '):
    self.neta=neta
    self.iettype=iettype
    self.ifltype=ifltype
    self.itetype=itetype
    self.isatype=isatype
    self.itrtype=itrtype


class ForcFreq(object):
  """
  Class of Boundary Forcing Frequency
  """
  def __init__(self,alpha=None,amig=None,ff=None,face=None):
    self.alpha=alpha
    self.amig=amig
    self.ff=ff
    self.face=face

	
class ConstPotentiel(object):
  """
  Class of Tidal Constituant Potential
  """
  def __init__(self,talpha=None,jspc=None,tamp=None,tfreq=None,tnf=None,tear=None):
    self.talpha=talpha
    self.jspc=jspc
    self.tamp=tamp
    self.tfreq=tfreq
    self.tnf=tnf
    self.tear=tear

def suppr_ligne_vide(fname):
    import os
    f1=open(fname)
    f2=open('bctides_tmp','wb')
    lines=f1.readlines()
    for line in lines:
        if line!='\r\n' and line!='\n':
            f2.write(str(line))
    f2.close()
    os.system('mv bctides_tmp '+str(fname)) 
    #print 'Empty lines suppressed '  


def read_bctides(bc_file):
    """
    Read a bctides.in file and return all the fields of a bctides including:
        date
        ntip
        list of tidal constituant potential
        nbfr
        list of boundary forcing frequencies
        list of open boundaries flags
    """
    suppr_ligne_vide(bc_file)
    f=open(bc_file)
    header=f.readline()
    bc=BcTides(date=header) # date bctide

    ntip, tip_dp=f.readline().split()[0:2]
    bc.tip_dp=tip_dp # attributing cut-off depth
    list_ntip=[]  # potential contituants
    for i in range(int(ntip)):
        talpha = f.readline().split()[0]
        jspc, tamp, tfreq, tnf, tear = f.readline().split()[:5]
        list_ntip.append(ConstPotentiel(talpha,jspc,tamp,tfreq,tnf,tear))
    bc.list_tip=list_ntip # attributing tip list

    nbfr=f.readline().split()[0]
    list_nbfr=[]  # forcing frequencies
    for i in range(int(nbfr)):
        alpha = f.readline().split()[0]
        amig, ff, face = f.readline().split()[:3]
        list_nbfr.append(ForcFreq(alpha,amig,ff,face))
    bc.list_bfr=list_nbfr # attributing bfr list

    # reading number of open boudaries
    nop=f.readline().split()[0]
    bc.nop=nop

    # reading list of open boundaries
    bc.list_openBound=[]
    # reading flags
    
    for i in range(int(nop)):	
        flag_l=f.readline().split()
        print(len(flag_l))
#		if len(flag_l)==5:
#			neta,iettype,ifltype,itetype,isatype=flag_l
#			flags=BcFlag(neta,iettype,ifltype,itetype,isatype)
#		if len(flag_l)==6:
#			neta,iettype,ifltype,itetype,istype,itrtype=flag_l
#			flags=BcFlag(neta,iettype,ifltype,itetype,isatype,itrtype)
#		bc.flags=flags
#		bc.list_openBound.append(flags)	
#		
    return bc


def read_bound_gr3(file):
    """
    Lecture des noeuds aux frontières ouvertes
    """
    print('read_bound_gr3... \n')
    f = open(file)
    fname = f.readline()
    print('hgrid file:', fname)
    nbelem, nbnode = f.readline().split()
    for i in range(int(nbelem) + int(nbnode)):
        f.readline()
    nbOpenBound = int(f.readline().split()[0])
    print('Number of open boundaries:', nbOpenBound)
    TotalNumb = int(f.readline().split()[0])
    print('Total number of open boundary nodes:', TotalNumb)
    list_bound = []
    Bound_num_list = arange(nbOpenBound) + 1
    count_bound = 0
    NumbBound = []
    for j in range(nbOpenBound):
        NumbBound.append(int(f.readline().split()[0]))
        print('Number of nodes for open boundary', Bound_num_list[count_bound], ':', NumbBound[j])
        for i in range(NumbBound[-1]):
            list_bound.append(int(f.readline().split()[0]))

        count_bound = count_bound + 1
        
    list_bound = array(list_bound)
    list_bound_coord = get_coord(list_bound, file)

    # print list_bound_coord
    return list_bound_coord, TotalNumb, nbOpenBound, NumbBound


def get_constant_name():
    """
    read contituant list name in FES2012 folder
    """
    list_name = []


def get_coord(list_bound, file):
    """
    Get coordinates for each node along the open open boundary
    Node numbers from read_bound_gr3
    """
    f = open(file)
    # list_coord=[]
    fname = f.readline()
    nbelem, nbnode = f.readline().split()
    indx = 1
    test_list = []
    list_coord = zeros([len(list_bound), 3])
    for i in range(int(nbnode)):
        numNode, x, y, z = f.readline().split()
        test = find(list_bound == int(numNode))
        if test or test == 0:
            test_list.append(test[0])  # ,int(numNode),float(x),float(y)])
            # list_coord.append([int(numNode),float(x),float(y)])
            list_coord[test] = [int(numNode), float(x), float(y)]
            # print i+2, indx, numNode, find(list_bound==int(numNode))
            indx = indx + 1
        # test_list=sort(test_list,axis=1)
        # list_coord.append([int(numNode),float(x),float(y)])

    # print test_list
    # list_coord=array(list_coord)
    return list_coord


def emprise(list_coord):
    """
    Définition de l'emprise nécessaire pour couvrir la grille gr3
    """
    up = round(list_coord[:, 2].max()) + 2
    down = round(list_coord[:, 2].min()) - 2
    # right=	round(list_coord[:,1].max())+362
    # left=	round(list_coord[:,1].min())+358
    right = round(list_coord[:, 1].max()) + 2
    left = round(list_coord[:, 1].min()) - 2
    limite = [up, down, left, right]
    return limite


def get_netcdf_coord(lat, lon, var):
    coord = []
    for i in range(len(lat)):
        for j in range(len(lon)):
            coord.append([lon[j], lat[i], var[i, j]])

    coord = array(coord)
    return coord


def read_constituant(const, limite):
    """
    Read Ha and Hg variables of constituant const
    """
    l = os.listdir(global_atlas_dir)
    for i in l:
        if i.split('_')[0].lower() == const.lower():
            const_name = i
    m = Dataset(path + const_name)
    print('Reading file ', const_name, ' ...')

    lat = m.variables['lat'][:]
    lon = m.variables['lon'][:]

    ind_lat = [find(lat > limite[1])[0], find(lat < limite[0])[-1]]
    ind_lon = [find(lon > limite[2])[0], find(lon < limite[3])[-1]]

    lat = lat[ind_lat[0]:ind_lat[1]]
    lon = lon[ind_lon[0]:ind_lon[1]]

    Ha = m.variables['Ha'][ind_lat[0]:ind_lat[1], ind_lon[0]:ind_lon[1]]
    Hg = m.variables['Hg'][ind_lat[0]:ind_lat[1], ind_lon[0]:ind_lon[1]]
    Ha = get_netcdf_coord(lat, lon, Ha)
    Hg = get_netcdf_coord(lat, lon, Hg)
    return Ha, Hg


def interpolate2gr3(const, list_bound_coord):
    """
    Interpolation of amplitude (Ha) and Phaselag (Hg) along boundaries
    nodes given in bound_coord list read in /Users/gael/Documents/Campagnes/FES2012/:
    """
    from scipy.interpolate import griddata
    limite = emprise(list_bound_coord)
    Ha, Hg = read_constituant(const, limite)
    # delete nan values in FES2012 data to avoid nan in interpolation
    Ha = delete(Ha, find(isnan(Ha[:, 2])), axis=0)
    Hg = delete(Hg, find(isnan(Hg[:, 2])), axis=0)
    # Tidal elevation amplitude
    Ha_bound = griddata((Ha[:, 0], Ha[:, 1]), Ha[:, 2], (list_bound_coord[:, 1], list_bound_coord[:, 2]),
                        method='nearest')
    if find(isnan(Ha_bound) == True).any():
        print('\t nan detected in Ha_bound')

    # Tidal elevation phaselag
    Hg_bound = griddata((Hg[:, 0], Hg[:, 1]), Hg[:, 2], (list_bound_coord[:, 1], list_bound_coord[:, 2]),
                        method='nearest')
    if find(isnan(Ha_bound) == True).any():
        print('\t nan detected in Hg_bound')

    # convertint amplitude (cm) into (m)
    Ha_bound = Ha_bound / 100.

    # phase positive
    Hg_bound = np.where(Hg_bound < 0, Hg_bound + 360, Hg_bound)

    return Ha_bound, Hg_bound


def WGS2UTM(x1, y1):
    inProj = Proj(init='epsg:4326')  # WGS84
    outProj = Proj(init='epsg:32620')  # UTM20

    x2, y2 = transform(inProj, outProj, x1, y1)
    coord = transpose([x2, y2])

    return coord


def UTM2WGS(x1, y1):
    inProj = Proj(init='epsg:32620')  # UTM20
    outProj = Proj(init='epsg:4326')  # WGS84

    x2, y2 = transform(inProj, outProj, x1, y1)

    return x2, y2


def write_bctides(gr3_file, bctide_file, out_file):
    """
    Writing a bctide file from a bctide class
    """
    class_bctides = read_bctides(bctide_file)

    list_coord, TotalNeta, nop, NumbBound = read_bound_gr3(gr3_file)
    # modifying number of open boundaries
    class_bctides.nop = nop

    # adapt number of open boundaries
    # default flags are elevation forcing at boundaries ( 3 0 0 0 )
    class_bctides.list_openBound = []
    for i in range(nop):
        class_bctides.list_openBound.append(BcFlag(NumbBound[i], 3, 0, 0, 0))

    print('Modifying number of node along open boundaries:', TotalNeta)

    f = open(out_file, 'w')			# Will destroy previous file.
    f.write(str(class_bctides.date))  	# writing date

    # writing tidal potential constituants
    f.write(str(len(class_bctides.list_tip)) + ' ' + class_bctides.tip_dp + ' !ntip, tip_dp \n')  # ntip and tip_dp
    for i in class_bctides.list_tip:
        f.write('\t' + i.talpha)
        f.write('\n')
        f.write(str(i.jspc) + ' ' + str(i.tamp) + ' ' + str(i.tfreq) + ' ' + str(i.tnf) + ' ' + str(i.tear) + '\n')

    # writing boundary forcing frequencies
    f.write(str(len(class_bctides.list_bfr)) + '  !nbfr \n')  # nbfr
    for i in class_bctides.list_bfr:
        f.write('\t' + i.alpha)
        f.write('\n')
        f.write(str(i.amig) + ' ' + str(i.ff) + ' ' + str(i.face) + '\n')

    f.write(str(class_bctides.nop) + ' !Number of Open Boundaries\n')

    # writing nb of boundary node and flags

    deb = 0
    fin = 0
    for i in range(nop):
        class_bctides.list_openBound[i].neta = NumbBound[i]
        fin = fin + int(class_bctides.list_openBound[i].neta)

        f.write(str(class_bctides.list_openBound[i].neta) + ' ' + str(class_bctides.list_openBound[i].iettype) \
                + ' ' + str(class_bctides.list_openBound[i].ifltype) + ' ' + str(
            class_bctides.list_openBound[i].itetype) \
                + ' ' + str(class_bctides.list_openBound[i].isatype) + ' ' + str(
            class_bctides.list_openBound[i].itrtype) + '\n')

        print('\nComputing open boundarie number:', i + 1)

        # writing amplitude and phase at boundary node for each constituant for each boundary segment
        for k in range(len(class_bctides.list_bfr)):
            print('Computing constituents ', class_bctides.list_bfr[k].alpha)
            Ha, Hb = interpolate2gr3(class_bctides.list_bfr[k].alpha, list_coord[deb:fin])
            f.write(str(class_bctides.list_bfr[k].alpha) + '\n')
            for j in range(int(class_bctides.list_openBound[i].neta)):
                f.write(str('%.12f' % Ha[j]) + '\t' + str('%.12f' % Hb[j]) + '\n')
        deb = deb + int(class_bctides.list_openBound[i].neta)

    f.close()


def suppr_ligne_vide(fname):
    import os
    f1 = open(fname)
    f2 = open('bctides_tmp', 'wb')
    lines = f1.readlines()
    for line in lines:
        if line != '\r\n' and line != '\n':
            f2.write(str(line))
    f2.close()
    os.system('mv bctides_tmp ' + str(fname))


# print 'Empty lines suppressed'
#################################################################################

write_bctides(hgrid_file_name, bctide_example, bctide_outfile)
