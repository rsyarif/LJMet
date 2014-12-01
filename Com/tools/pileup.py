#!/usr/bin/env python
#########################################################################
#
# pileup.py
#
# Module for pileup-related procedures
# - save MC pileup scenario from CMSSW to a histogram in a file
# - prepare all needed for pileup reweighting for a list of datasets
# 
#
# Usage:
#        Grab a PU scenario from CMSSW:
#        ./analyze.py --mc-pileup <scenario>
#
#
# Authors:
#          Gena Kukartsev, July 2012
#
#
#########################################################################



import os
import glob
import re
import sys
import optparse
import ConfigParser
import json_manip
import subprocess
import datetime

from inspect import getmembers



legend = '[pileup]:'



def GetPuFiles(package):
    
    #get the release base location
    cmssw_release_base = os.environ['CMSSW_RELEASE_BASE']

    # path
    _path = cmssw_release_base+'/src/'+package+'/python'

    # get lists of py files
    _files = glob.glob(_path+'/*.py')

    return _files



def GetPuScenarios(files):

    _scenarios = []
    for path in files:
        _scenario = os.path.basename(path)
        _scenario = re.sub('_cff.py', '', _scenario)
        _scenario = re.sub('_cfi.py', '', _scenario)
        _scenario = re.sub('.py', '', _scenario)
        if _scenario.find('_') <= 0:
            continue
        # remove up to the first underscore
        #_scenario = _scenario[_scenario.find('_')+1:]
        _scenarios.append(_scenario)
        #print legend, _scenario

    return _scenarios



def GetMatchingScenarios(scenarios, pattern):
    #
    # find all PU scenarios that match pattern
    #
    _match = []
    for sc in scenarios:
        if sc.find(pattern)>=0:
            _match.append(sc)

    return _match



def RecursiveHasAttr(obj, attr):

    try:
        left, right = attr.split('.', 1)
    except:
        return hasattr(obj, attr)

    return RecursiveHasAttr(getattr(obj, left, None), right)



def RecursiveGetAttr(obj, attr, default=None):
    
    try:
        left, right = attr.split('.', 1)
    except:
        return getattr(obj, attr, default)
    
    return RecursiveGetAttr(getattr(obj, left, default), right, default)



def SavePuHist(files, sc, full_sim, prefix = ''):
    #
    # saves the 'sc' scenario to a histogram
    #

    # possible names of the objects which contain PU info
    mod_names = []
    mod_names.append('famosPileUp.PileUpSimulator')
    mod_names.append('mix.input.nbPileupEvents')
    mod_names.append('mixGenPU.input.nbPileupEvents')

    _filename = 'pileup_'
    if full_sim:
        _filename += 'FULLSIM_'
    else:
        _filename += 'FASTSIM_'
    _filename = prefix + _filename + sc+'.root'


    _nfound = 0
    for file in files:
        if file.find(sc+'_cfi')>=0:
            _nfound += 1
            _module = re.sub('.py', '', os.path.basename(file))
            sys.path.append(os.path.dirname(file))

            exec('import '+_module+' as pu_mod')

            for mname in mod_names:
                if RecursiveHasAttr(pu_mod, mname):
                    _mod = RecursiveGetAttr(pu_mod, mname)

            _bin  = _mod.probFunctionVariable
            _prob = _mod.probValue

            _nbins = len(_bin)
            import ROOT
            _file = ROOT.TFile(_filename, 'recreate')
            _hist = ROOT.TH1F('MC_pileup', 'MC_pileup', _nbins, 0, _nbins)
            for bin in _bin:
                _hist.SetBinContent(bin+1, _prob[bin])

            print legend, 'saving to file', _filename
            _file.Write()

    if _nfound < 1:
        print legend, 'WARNING: PileUp information not found, cannot create a histogram!'

    return



    
def GetCvsHistogram(options, prefix = ''):
    
    fastsim_cvs = 'FastSimulation/PileUpProducer'
    fullsim_cvs = 'SimGeneral/MixingModule'

    fastsim_files = GetPuFiles(fastsim_cvs)
    fullsim_files = GetPuFiles(fullsim_cvs)


    fastsim_scenarios = GetPuScenarios(fastsim_files)
    fullsim_scenarios = GetPuScenarios(fullsim_files)

    if options.mc_pileup != 'ALL':

        if options.fastsim_exact:
            fastsim_scenarios = [options.mc_pileup]
            fullsim_scenarios = []

        elif options.fullsim_exact:
            fullsim_scenarios = [options.mc_pileup]
            fastsim_scenarios = []

        else:
            fastsim_scenarios = GetMatchingScenarios(fastsim_scenarios,
                                                     options.mc_pileup)
            fullsim_scenarios = GetMatchingScenarios(fullsim_scenarios,
                                                     options.mc_pileup)


    if options.mc_pileup == 'FILES':
        for _file in fastsim_files:
            print legend, 'FASTSIM:', _file
        for _file in fullsim_files:
            print legend, 'FULLSIM:', _file


    if ( len(fastsim_scenarios) + len(fullsim_scenarios) ) < 1:
        print legend, 'no matching pileup scenarios found'
        return

    elif ( len(fastsim_scenarios) + len(fullsim_scenarios) ) > 1:
        print legend, 'more than one PU scenario found'
        print legend, 'possible FastSim pileup scenarios:'
        for sc in fastsim_scenarios:
            print legend, sc
        print legend, 'possible FullSim pileup scenarios:'
        for sc in fullsim_scenarios:
            print legend, sc

        return


    if len(fastsim_scenarios)==1:
        print legend, 'one FastSim pileup scenario found'
        SavePuHist(fastsim_files, fastsim_scenarios[0], full_sim = False,
                   prefix = prefix)

    if len(fullsim_scenarios)==1:
        print legend, 'one FullSim pileup scenario found'
        SavePuHist(fullsim_files, fullsim_scenarios[0], full_sim = True,
                   prefix = prefix)

    return



class PileupManager:

    def __init__(self, fname):
        
        print legend, 'initializing Pileup manager'

        self.config = None
        self.conf_map = {} # map to store config conveniently for random access

        # read config file
        with open(fname) as file:
            
            self.config = ConfigParser.ConfigParser()
            self.config.readfp(file)


            # populate self.conf_map
            for section in self.config.sections():

                if section not in self.conf_map:
                    self.conf_map[section] = {}

                for item in self.config.items(section):
                    self.conf_map[section][item[0]] = item[1]

        # get CMSSW base dir
        self.cmssw_base = os.environ['CMSSW_BASE']
        self.scram_arch = os.environ['SCRAM_ARCH']
        print legend, 'release base dir:', self.cmssw_base

        # package name
        self.package_name = 'LJMet/Com'

        # file locations
        self.json_location   = 'data/json'
        self.pileup_location = 'data/pileup'

        # timestamp
        self.now = datetime.datetime.now()
        self.timestamp = self.now.strftime('_%d%b%Y').lower()

        # provenance from config file
        # check if config info is available
        if not self.config:
            print legend, 'config file is not parsed yet'
            sys.exit(-1)

        # check if PROVENANCE section exists
        if 'PROVENANCE' not in self.conf_map:
            print legend, 'No PROVENANCE section found in config file'
            sys.exit(-1)

        # check if name is defined
        #_items = self.config.items('PROVENANCE')
        #print self.conf_map['PROVENANCE']
        if 'prefix' not in self.conf_map['PROVENANCE']:
            print legend, 'No prefix defined in config file'
            sys.exit(-1)

        # check if version is defined
        if 'version' not in self.conf_map['PROVENANCE']:
            print legend, 'No version defined in config file'
            sys.exit(-1)

        self.prefix  = self.conf_map['PROVENANCE']['prefix']
        self.version = self.conf_map['PROVENANCE']['version']

        



    def GetLegend(self):

        _legend  = '['+legend.strip().strip(':').strip('[]')+'.'+self.__class__.__name__+']:'

        return _legend



    def CheckFileExists(self,fname):
        #
        # check if file exists in the release
        #

        #_path = self.cmssw_base+'/src/'+self.package_name+'/'+fname.strip().strip('/')

        return os.path.exists(fname)


        
    def GetFullPath(self,fname):
        #
        # add base and package names to path
        #

        return self.cmssw_base+'/src/'+self.package_name+'/'+fname.strip().strip('/')

        
    def ReadConfigFile(self,fname):
        
        with open(fname) as file:
            
            self.config = ConfigParser.ConfigParser()
            self.config.readfp(file)


            # populate self.conf_map
            for section in self.config.sections():

                if section not in self.conf_map:
                    self.conf_map[section] = {}

                for item in self.config.items(section):
                    self.conf_map[section][item[0]] = item[1]



    def GetMcPileupScenario(self):

        # check if config info is available
        if not self.config:
            print legend, 'config file is not parsed yet'
            return -1

        # check if PILEUP section exists
        if 'PILEUP' not in self.conf_map:
            print legend, 'No PILEUP section found in config file'
            return -1

        # check if scenario is defined
        _items = self.config.items('PILEUP')
        if 'mc_scenario' not in self.conf_map['PILEUP']:
            print legend, 'No MC scenario defined in config file'
            return -1

        _scenario   = self.conf_map['PILEUP']['mc_scenario'].split(';')[0].strip()
        _is_fullsim = self.conf_map['PILEUP']['mc_scenario'].split(';')[1].strip()=='fullsim'


        # configure options
        _prefix = 'pileup_'
        _options = optparse.Values()
        _options.ensure_value('mc_pileup', _scenario)
        if _is_fullsim:
            _options.ensure_value('fullsim_exact', True)
            _options.ensure_value('fastsim_exact', False)
            _prefix += 'FULLSIM_'
        else:
            _options.ensure_value('fullsim_exact', False)
            _options.ensure_value('fastsim_exact', True)
            _prefix += 'FASTSIM_'


        #print _options.mc_pileup

        
        # check if PU hist file exists
        _baseprefix = self.cmssw_base+'/src/'+self.package_name+'/'+self.pileup_location+'/'
        _filename = _baseprefix+_prefix+_scenario+'.root'
        if self.CheckFileExists(_filename):
            print legend, _filename, 'already exists'
            return
        else:
            GetCvsHistogram(_options, prefix = _baseprefix)
            print legend, _filename, 'created'

        return
            


    def GetJsonFromDbs(self, dataset, dsvar=None):
        #
        # Query DBS for the dataset,
        # Make JSON file
        #
        # dsvar contains dataset var name in the config file
        # we get non-default DBS instance from there
        #

        legend = self.GetLegend()

        print legend, 'creating JSON for', dataset

        json_name = self.GetFullPath(self.json_location+'/'
                                     +dataset.strip('/').replace('/', '_').replace('-', '_')
                                     +'_DBS_'+self.prefix
                                     +'_v'+self.version+'.json')

        
        # check if already exists
        #_fullpath = self.GetFullPath(json_name)
        if self.CheckFileExists(json_name):
            print self.GetLegend(), json_name, 'already exists'
            return json_name


        # create new JSON
        print self.GetLegend(), json_name, 'being created...'
        _options = optparse.Values()
        _dataset = []
        _dataset.append(dataset)
        _in_file = []
        _options.ensure_value('dataset', _dataset)
        _options.ensure_value('in_file', _in_file)
        _options.ensure_value('out_file', json_name)
        _options.ensure_value('debug', True)
        _options.ensure_value('intersect', False)

        _dbs_instance='http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'


        # get non-default DBS instance if specified
        if dsvar:
            if len(dsvar.split(';'))>1:
                _dbs_alias = dsvar.split(';')[1].strip()
                if _dbs_alias=='phys02':
                    _dbs_instance='http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'
                    #print 'DEBUG!!!', _dbs_instance

        _options.ensure_value('dbs_instance', _dbs_instance)
        #print _options

        json_manip.GetJsonFile(_options)
        

        return json_name



    def MakeCustomJson(self, dataset, dsvar=None):
        #
        # Query DBS for the dataset,
        # compare to official JSON,
        # Make JSON that is intersection with official JSON
        #
        # dsvar contains name of dataset in config file
        #

        legend = self.GetLegend()

        
        # get JSON from DBS
        _json_dbs = self.GetJsonFromDbs(dataset,dsvar)


        json_name = self.GetFullPath(self.json_location+'/'
                                     +dataset.strip('/').replace('/', '_').replace('-', '_')
                                     +'_'+self.prefix
                                     +'_v'+self.version+'.json')

        #print self.GetLegend(), json_name

        
        # check if already exists
        #_fullpath = self.GetFullPath(json_name)
        if self.CheckFileExists(json_name):
            print self.GetLegend(), json_name, 'already exists'
            return json_name


        # check if JSON section exists
        if 'JSON' not in self.conf_map:
            print legend, 'No JSON section found in config file'
            return -1

        # check if official json is specified
        if 'official_json' not in self.conf_map['JSON']:
            print legend, 'Official JSON file not specified in config file'
            return -1
        
        # check if official JSON exists
        _official_json = self.GetFullPath(self.json_location+'/'+self.conf_map['JSON']['official_json'])
        if not self.CheckFileExists(_official_json):
            print legend, 'Official JSON file not found!'
            return -1

        # create new JSON
        print self.GetLegend(), json_name, 'being created...'
        _options = optparse.Values()
        _dataset = []

        #debug
        #_dataset.append(dataset)
        _in_file = []
        
        #debug
        _in_file.append(_json_dbs)

        _in_file.append(_official_json)
        _options.ensure_value('dataset', _dataset)
        _options.ensure_value('in_file', _in_file)
        _options.ensure_value('out_file', json_name)
        _options.ensure_value('debug', True)
        _options.ensure_value('intersect', True)

        #debug
        #_dbs_instance='http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'
        #_dbs_instance='http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
        #_options.ensure_value('dbs_instance', _dbs_instance)
        #print _options

        json_manip.GetJsonFile(_options)
        

        return json_name



    def MakePileupHist(self, json, outfile):
        #
        # Run PileupCalc.py and generate pileup
        # histogram for a JSON file
        #

        

        # check if JSON section exists
        if 'JSON' not in self.conf_map:
            print legend, 'No JSON section found in config file'
            return -1

        # check if lumi json is specified
        if 'lumi_json' not in self.conf_map['JSON']:
            print legend, 'Lumi JSON file not specified in config file'
            return -1
        
        # check if lumi JSON exists
        _lumi_json = self.GetFullPath(self.json_location+'/'+self.conf_map['JSON']['lumi_json'])
        if not self.CheckFileExists(_lumi_json):
            print legend, 'Lumi JSON file not found!'
            return -1


        _min_bias_xsec = None
        # check if PILEUP section exists
        if 'PILEUP' not in self.conf_map:
            print legend, 'No PILEUP section found in config file'
        else:
            if 'min_bias_xsec' not in self.conf_map['PILEUP']:
                print legend, 'No Minbias cross section specified in config file'
            else:
                _min_bias_xsec = self.conf_map['PILEUP']['min_bias_xsec'].strip()

        if _min_bias_xsec:
            print legend, 'using minbias cross section:', _min_bias_xsec
        else:
            _min_bias_xsec = '69300'
            print legend, 'minbias cross section not specified, using', _min_bias_xsec


        # run PileupCalc.py
        _pileup_calc = self.cmssw_base+'/bin/'+self.scram_arch+'/'+'pileupCalc.py'
        _pipe = subprocess.Popen(['python', _pileup_calc,
                                  '-i', json,
                                  '--inputLumiJSON', _lumi_json,
                                  '--calcMode', 'true',
                                  #'--minBiasXsec', '73500',
                                  '--minBiasXsec', _min_bias_xsec,
                                  '--maxPileupBin', '60',
                                  '--numPileupBins', '60',
                                  outfile],
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
        
        stdout, stderr = _pipe.communicate()


        return



    def GetLumi(self, json):
        #
        # Get lumi overview from lumiCalc2.py
        #

        # run PileupCalc.py
        _lumi_calc = self.cmssw_base+'/bin/'+self.scram_arch+'/'+'lumiCalc2.py'
        _pipe = subprocess.Popen(['python', _lumi_calc,
                                  '-b', 'stable',
                                  '-i', json,
                                  'overview'],
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
        
        stdout, stderr = _pipe.communicate()

        return stdout

        

    def GetDataPileup(self):
        #
        # Get pileup histogram for given datasets
        # - get JSON for each dataset
        # - merge data JSON
        # - get PU
        #

        # for each dataset, create JSON and pileup hist
        _json_list = [] # list of all json to merge later
        for ds in self.conf_map['DATA']:
            _ds_name = self.conf_map['DATA'][ds]

            # debug
            #print ds, _ds_name

            # create dataset JSON file
            #_json_dbs = self.GetJsonFromDbs(_ds_name)
            _json = self.MakeCustomJson(_ds_name, ds)
            _json_list.append(_json)

            # create dataset pileup hist
            _pu_name = self.GetFullPath(self.pileup_location+'/'
                                        +'pileup_'
                                        +_ds_name.strip('/').replace('/', '_').replace('-', '_')
                                        +'_'+self.prefix
                                        +'_v'+self.version+'.root')


            print self.GetLegend(), 'creating pileup histogram for', ds, '...'
            if not self.CheckFileExists(_pu_name):
                self.MakePileupHist(_json, _pu_name)

            # create dataset lumi overview
            _lumi_name = self.GetFullPath(self.json_location+'/'
                                          +_ds_name.strip('/').replace('/', '_').replace('-', '_')
                                          +'_'+self.prefix
                                          +'_v'+self.version+'.lumi')
            print self.GetLegend(), 'getting luminosity info for', ds, '...'
            if not self.CheckFileExists(_lumi_name):
                _stdout = self.GetLumi(_json)
                with open(_lumi_name, 'w') as _lumi_file:
                    _lumi_file.write(_stdout)
                    print self.GetLegend(), _lumi_name, 'created'
            

        # merge all dataset JSONs into master JSON
        #_json_master = self.cmssw_base+'/src/'+self.package_name+'/'+self.json_location+'/'+'master'+self.timestamp+'.json'
        _json_master = self.cmssw_base+'/src/'+self.package_name+'/'+self.json_location+'/'+self.prefix+'_v'+self.version+'.json'
        print self.GetLegend(), 'JSON master file being created...'
        _options = optparse.Values()
        _in_file = _json_list
        _dataset=[]
        _options.ensure_value('in_file', _in_file)
        _options.ensure_value('dataset', _dataset)
        _options.ensure_value('out_file', _json_master)
        _options.ensure_value('debug', True)
        _options.ensure_value('intersect', False)
        #print _options
        json_manip.GetJsonFile(_options)


        # make master data pileup histogram
        print self.GetLegend(), 'pileup master file being created...'
        _pu_master = self.cmssw_base+'/src/'+self.package_name+'/'+self.pileup_location+'/'+'pileup_'+self.prefix+'_v'+self.version+'.root'
        if not self.CheckFileExists(_pu_master):
            self.MakePileupHist(_json_master, _pu_master)
            print self.GetLegend(), _pu_master, 'created'
        else:
            print self.GetLegend(), _pu_master, 'already exists'


        # get overall lumi overview
        _lumi_master = self.cmssw_base+'/src/'+self.package_name+'/'+self.json_location+'/'+self.prefix+'_v'+self.version+'.lumi'
        print self.GetLegend(), 'getting luminosity info for all data...'
        if not self.CheckFileExists(_lumi_master):
            _stdout = self.GetLumi(_json_master)
            with open(_lumi_master, 'w') as _lumi_file:
                _lumi_file.write(_stdout)
                print self.GetLegend(), _lumi_master, 'created'

        return



    def RunPileupJob(self):

        # get MC pileup hist
        self.GetMcPileupScenario()

        #self.MakeCustomJson(self.conf_map['DATA']['ds_run2012c_part1'])
        self.GetDataPileup()



def RunPileupJob(options):

    if not options.config_file:

        print legend, 'config file not specified, exiting'
        sys.exit(-1)


    pu = PileupManager(options.config_file)
    pu.RunPileupJob()
        



#########################################################################
#
# Main
#


if __name__=='__main__':

    legend = '[pileup]:'



########################################
#
# banner
#
    def banner():
        print '''
+--------------------------------------------------------------
|
| pileup.py
|
| Pileup info and more
|
| authors:
|           Gena Kukartsev, July 2012
|
+--------------------------------------------------------------
    '''
    banner()



    import sys



########################################
#
# parse command line parameters
#
    from optparse import OptionParser
    add_help_option = '''
./pileup.py -c config.cfg
'''

    parser = OptionParser(add_help_option)


    parser.add_option("-c", "--config-file", dest="config_file", default=None,
                  help="Name of the config file")

    
    print legend, 'parsing command line options...',
    (options, args) = parser.parse_args()
    print 'done'


    # run
    RunPileupJob(options)

