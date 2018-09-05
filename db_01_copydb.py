import pdb
import os,glob,sys
group = sys.argv[1]
path = '/data/shangzhong/Proteogenomics/ppl'
msgf = '/home/shangzhong/Codes/Proteogenomics/MSGF_Enosi/MSGFPlus.jar'
# path = '/media/lewislab/Dropbox (UCSD SBRG)/users/shangzhong/Pipeline'

group = str(group)
db_path = path + '/DBlist'
spectrum_path = path + '/group' + str(group) + '/convert_spectrum'
spectrum_files = sorted(glob.glob(spectrum_path + '/*.mgf'))
list_path = path + '/group' + group +  '/list'
if not os.path.exists(list_path): os.mkdir(list_path)

if group == '1':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '2':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '3':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '4':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '5':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '6':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '7':
	params = '-t 15ppm -tda 1 -m 3 -inst 0 -e 1 -ti 0,1 -ntt 2'
elif group == '8':
	params = '-t 20ppm -tda 1 -m 3 -inst 3 -e 1 -ti 0,1 -ntt 2'
elif group == '9':
	params = '-t 10ppm -tda 1 -m 0 -inst 3 -e 1 -ti 0,1 -ntt 1'


all_dbs = []
dbs = sorted(glob.glob(db_path+'/*'))
for db in dbs:
	sub_path = db + '/prep'
	sub_dbs = glob.glob(sub_path + '/*')
	for sub_db in sub_dbs:		
		all_dbs.append(glob.glob(sub_db+'/*.fa')[0])

cmds = []
for mgf in spectrum_files:
	mgf = '/'.join(mgf.split('/')[-2:])
	for db in all_dbs:
		sequence = db.split('/')[-4]
		fn = mgf.split('/')[-1][:-4]+'_'+db.split('/')[-1][:-3]
		out_fn='MSGFresult/{seq}/{fn}.mzid'.format(seq=sequence,fn=fn)
		if os.path.exists(otu_fn): continue
		cmd = ('java -Xmx10G -jar {msgf} '
			'-s {mgf} -d \'{db}\' -o MSGFresult/{seq}/{fn}.mzid {param} '
			'-mod ../mod/group{g}.txt -thread 1').format(mgf=mgf,db=db,param=params,seq=sequence,fn=fn,
			g=str(group),msgf=msgf)
		cmds.append(cmd)
with open(list_path+'/LaunchMSGF.txt','w') as out:
	out.write('\n'.join(cmds) + '\n')