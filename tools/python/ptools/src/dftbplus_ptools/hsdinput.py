#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""
Module for interaction with the hsd input file.
"""


import os
import hsd
import numpy as np

class Changehsd:
    """Class for changing the hsd input"""

    def __init__(self, filename="dftb_in.hsd", directory=".", dictionary=None):
        """Initialises the changehsd class

        Args:
            directory (str): directory of file
            filename (str): name of file
            dictionary (dict): option to change existing dict
        """
        self._filename = filename
        self._directory = directory
        if dictionary is None:
            self._dictionary = self._read_hsd()
        else:
            self._dictionary = dictionary


    def _read_hsd(self):
        """Reads dictionary from file

        Returns:
            (dict): dictionary contained in file, empty if file not found/
                readable
        """
        path = os.path.join(self._directory, self._filename)
        try:
            return hsd.load(path)
        except OSError:
            return {}


    def set_hsd(self, dictionary):
        """Sets self._dictionary

        Args:
            dictionary (dict): dictionary to change
        """
        self._dictionary = dictionary


    def write_hsd(self):
        """Writes dictionary to file

        Returns:
            (file): file named 'filename' at 'directory'
        """
        path = os.path.join(self._directory, self._filename)
        hsd.dump(self._dictionary, path)


    def get_hsd(self):
        """Function for returning self._dictionary

        Returns:
            (dict): contains hsd dictionary
        """
        return self._dictionary


    def write_resultstag(self, value=True):
        """Function for setting 'WriteResultsTag'

        Args:
            value (bool): True for writing results.tag
        """
        self._set_keyword('Options')
        self._dictionary['Options']['WriteResultsTag'] = value


    def calc_forces(self, value=True):
        """Function for setting 'CalculateForces'

        Args:
            value (bool): True for calculating forces
        """
        self._set_keyword('Analysis')
        self._dictionary['Analysis']['CalculateForces'] = value


    def calc_charges(self, value=True):
        """Function for setting 'MullikenAnalysis'

        Args:
            value (bool): True for calculating Charges
        """
        self._set_keyword('Analysis')
        self._dictionary['Analysis']['MullikenAnalysis'] = value


    def set_driver(self, driver, drivermaxforce=None, drivermaxsteps=None):
        """Function for setting 'Driver', 'MaxForceComponent' and 'MaxSteps'

        Args:
            driver (str): DFTB+ driver for geometry optimization
            drivermaxforce (float): max. force component as convergence
                criterion of geometry optimization
            drivermaxsteps (int): max. number of geometry steps
        """
        self._dictionary['Driver'] = {}
        self._dictionary['Driver'][str(driver)] = {}
        if drivermaxforce is not None:
            self._dictionary['Driver'][str(driver)]['MaxForceComponent'] = \
                drivermaxforce
        if drivermaxsteps is not None:
            self._dictionary['Driver'][str(driver)]['MaxSteps'] = \
                drivermaxsteps


    def set_scc(self, value=True):
        """Function for setting 'self-consistent calculations'

        Args:
            value (bool): True for self-consistent calculations
        """
        self._set_keyword('Hamiltonian', 'DFTB')
        self._dictionary['Hamiltonian']['DFTB']['Scc'] = value


    def set_scctol(self, value=1E-005):
        """Function for setting 'convergence criterion of SCC cycles'

        Args:
            value (float): convergence criterion of SCC cycles
        """
        self._set_keyword('Hamiltonian', 'DFTB')
        self._dictionary['Hamiltonian']['DFTB']['SccTolerance'] = value


    def set_skdir(self, skdir):
        """Function for setting 'Slater-Koster files'

        Args:
            skdir (str/dict): path(s) to Slater-Koster files

        Raises
            TypeError: if skdir has wrong Type
        """
        self._set_keyword('Hamiltonian', 'DFTB')
        if isinstance(skdir, str):
            if not skdir.endswith('/'):
                skdir += '/'
            skfiledict = {}
            skfiledict['Prefix'] = skdir
            skfiledict['Separator'] = '"-"'
            skfiledict['Suffix'] = '".skf"'

            self._dictionary['Hamiltonian']['DFTB']['SlaterKosterFiles'] = {}
            self._dictionary['Hamiltonian']['DFTB']['SlaterKosterFiles'] \
                ['Type2FileNames'] = skfiledict

        elif isinstance(skdir, dict):
            for key, path in skdir.items():
                skdir[key] = '"' + path.strip('"').strip("'") + '"'
            self._dictionary['Hamiltonian']['DFTB']['SlaterKosterFiles'] = \
                skdir
        else:
            raise TypeError('Unexpected object type: "skdir" ' +
                            f'is type "{type(skdir).__name__}"' +
                            ' instead of "dict" or "str"')


    def set_kpts(self, kpts):
        """Function for setting 'K-points'

        Args:
            kpts (list/tuple): K-points for Brillouin zone sampling
                (periodic structures)

        Raises:
            ValueError: if K-Point definition is illegal
        """
        # Handle different K-point formats
        # (note: the ability to handle bandpaths has not yet been
        # implemented)
        self._set_keyword('Hamiltonian', 'DFTB')
        if np.array(kpts).ndim == 1:
            # Case: K-points as (gamma-centered) Monkhorst-Pack grid
            mp_mesh = kpts
            offsets = [0.] * 3
            props = [elem for elem in mp_mesh if isinstance(elem, str)]
            props = [elem.lower() for elem in props]
            tgamma = 'gamma' in props and len(props) == 1
            if tgamma:
                mp_mesh = mp_mesh[:-1]
                eps = 1e-10
                for ii in range(3):
                    offsets[ii] *= mp_mesh[ii]
                    assert abs(offsets[ii]) < eps or \
                        abs(offsets[ii] - 0.5) < eps
                    if mp_mesh[ii] % 2 == 0:
                        offsets[ii] += 0.5
            elif not tgamma and len(mp_mesh) != 3:
                raise ValueError('Illegal K-Points definition: ' +
                                 str(kpts))
            kpts_mp = np.vstack((np.eye(3) * mp_mesh, offsets))
            self._dictionary['Hamiltonian']['DFTB']['KPointsAndWeights'] = {}
            self._dictionary['Hamiltonian']['DFTB']['KPointsAndWeights'] \
                ['SupercellFolding'] = kpts_mp

        elif np.array(kpts).ndim == 2:
            # Case: single K-points explicitly as (N x 3) or (N x 4) matrix
            kpts_coord = np.array(kpts)
            if np.shape(kpts_coord)[1] == 4:
                kptsweights = kpts_coord
                kpts_coord = kpts_coord[:, :-1]
            elif np.shape(kpts_coord)[1] == 3:
                kptsweights = np.hstack([kpts_coord, [[1.0], ] *
                                         np.shape(kpts_coord)[0]])
            else:
                raise ValueError('Illegal K-Points definition: ' +
                                 str(kpts))
            self._dictionary['Hamiltonian']['DFTB']['KPointsAndWeights'] = \
                kptsweights
        else:
            raise ValueError('Illegal K-Points definition: ' + str(kpts))


    def set_filling(self, filling, order=None, temperature=None):
        """Function for setting 'filling methode', 'filling order' and
        'temperature'

        Args:
            temperature (float): electron temperature in energy units
            filling (str): filling type for elecron levels
            filling_order (int): order of the Methessel-Paxton scheme

        Raises:
            ValueError: if filling methode unknown
        """
        self._set_keyword('Hamiltonian', 'DFTB')
        if filling.lower() == 'methfesselpaxton':
            self._dictionary['Hamiltonian']['DFTB']['Filling'] = {}
            self._dictionary['Hamiltonian']['DFTB']['Filling']\
                ['MethfesselPaxton'] = {}
            if temperature is not None:
                self._dictionary['Hamiltonian']['DFTB']['Filling']\
                    ['MethfesselPaxton']['Temperatur'] = temperature
            if order is not None:
                self._dictionary['Hamiltonian']['DFTB']['Filling']\
                    ['MethfesselPaxton']['Order'] = order
        elif filling.lower() == 'fermi':
            self._dictionary['Hamiltonian']['DFTB']['Filling'] = {}
            self._dictionary['Hamiltonian']['DFTB']['Filling']['Fermi'] = {}
            if temperature is not None:
                self._dictionary['Hamiltonian']['DFTB']['Filling']['Fermi']\
                    ['Temperatur'] = temperature
        else:
            raise ValueError(f'Unknown filling methode "{filling}"')


    def set_maxang(self, maxangs=None, try_reading=None):
        """Function for setting 'max. angular momenta'

        Args:
            maxang (dict/None): max. angular momentum of atom types
            try_reading (list): list of atom types with unknown max. angular
                momentum to try to read from Slater-Koster files

        Raises:
            ValueError: if Slater-Koster files not specified or if obtained
                max. angular momentum is out of range (s-f) or can't be read
                from Slater-Koster file or if element is in try_reading and
                maxangs
        """
        if maxangs is None:
            maxangs = {}

        if try_reading is not None:
            for element in try_reading:
                if element in maxangs:
                    raise ValueError(f"'{element}' is already specified in " +
                                     "'maxangs'!")


            try:
                self._dictionary['Hamiltonian']['DFTB']['SlaterKosterFiles']
            except KeyError:
                raise ValueError("Slater-Koster-Files not available. Please" +
                                 " set them using the function 'set_skdir'.")
            skdir = self._dictionary['Hamiltonian']['DFTB']\
                ['SlaterKosterFiles']
            try:
                prefix = skdir['Type2FileNames']['Prefix'].strip('"').strip(
                    "'")
            except KeyError:
                prefix = skdir

            for species in try_reading:
                if isinstance(prefix, str):
                    path = os.path.join(prefix, '{0}-{0}.skf'.format(species))
                else:
                    try:
                        path = skdir['{0}-{0}'.format(species)].replace('"',
                                                                        '')
                    except KeyError:
                        raise ValueError(f"No path for '{species}-{species}'" +
                                         " specified!")
                maxang = self.read_max_angular_momentum(path)

                if maxang is None:
                    msg = 'Error: Could not read max. angular momentum ' + \
                        'from Slater-Koster file ' + \
                        '"{0}-{0}.skf".'.format(species) + \
                        ' Please specify manually!'
                    raise ValueError(msg)

                if maxang not in range(0, 4):
                    msg = 'The obtained max. angular momentum from ' + \
                          'Slater-Koster file ' + \
                          '"{0}-{0}.skf"'.format(species) + \
                          ' is out of range. Please check!'
                    raise ValueError(msg)
                maxang = '"{}"'.format('spdf'[maxang])
                maxangs[species] = maxang

        self._set_keyword('Hamiltonian', 'DFTB')
        self._dictionary['Hamiltonian']['DFTB']['MaxAngularMomentum'] = maxangs


    def ignore_unprocessed_nodes(self, value=True):
        """Function for setting 'IgnoreUnprocessedNodes'

        Args:
            value (bool): True for ignoring unprocessed nodes
        """
        self._set_keyword('ParserOptions')
        self._dictionary['ParserOptions']['IgnoreUnprocessedNodes'] = value


    def set_geometry(self, geometry):
        """Function for setting 'geometry'

        Args:
            geometry (dftbplus_ptools.geometry.Geometry object): contains the
                geometry informations
        """
        specieslist = geometry.specienames
        indexes = geometry.indexes

        geodict = {}
        geodict['TypeNames'] = specieslist
        typesandcoords = np.empty([len(indexes), 4], dtype=object)
        typesandcoords[:, 0] = indexes
        typesandcoords[:, 1:] = np.array(geometry.coords, dtype=float)
        geodict['TypesAndCoordinates'] = typesandcoords
        geodict['TypesAndCoordinates.attrib'] = 'Angstrom'

        periodic = geometry.periodic

        if periodic:
            geodict['Periodic'] = 'Yes'
            geodict['LatticeVectors'] = geometry.latvecs
            geodict['LatticeVectors.attrib'] = 'Angstrom'
        else:
            geodict['Periodic'] = 'No'

        self._dictionary['Geometry'] = geodict


    def get_basic_input(self, skdir='./', driver=None, drivermaxforce=None,
                        drivermaxsteps=None, scc=False, scctol=None,
                        maxang=None, try_reading=None, kpts=None,
                        temperature=None, filling=None, filling_order=None,
                        resultstag=True, unprocessed=True, forces=True):
        """Generates a suitable Python dictionary for a calculation with the
        hamiltonian set to DFTB.

        Args:
            skdir (str/dict): path(s) to Slater-Koster files
            driver (str): DFTB+ driver for geometry optimization
            drivermaxforce (float): max. force component as convergence
                criterion of geometry optimization
            drivermaxsteps (int): max. number of geometry steps
            scc (bool): True for self-consistent calculations
            scctol (float): convergence criterion of SCC cycles
            maxang (dict): max. angular momentum of atom types
            try_reading (list): list of atom types with unknown max. angular
                momentum to try to read from Slater-Koster files
            kpts (list/tuple): K-points for Brillouin zone sampling
                (periodic structures)
            temperature (float): electron temperature in energy units
            filling (str): filling type for elecron levels
            filling_order (int): order of the Methessel-Paxton scheme
            resultstag (bool): True for writing results.tag
            unprocessd (bool): True for ignoring unprocessed nodes
            forces (bool): True for calculating forces
        """

        self.set_skdir(skdir)
        self.set_driver(driver, drivermaxforce, drivermaxsteps)
        self.set_scc(scc)
        if scctol is not None:
            self.set_scctol(scctol)
        if maxang is not None:
            self.set_maxang(maxang, try_reading)
        else:
            self.set_maxang(try_reading=try_reading)
        if kpts is not None:
            self.set_kpts(kpts)
        if filling is not None:
            self.set_filling(filling, filling_order, temperature)
        self.write_resultstag(resultstag)
        self.ignore_unprocessed_nodes(unprocessed)
        self.calc_forces(forces)


    def _set_keyword(self, keyword1, keyword2=None, keyword3=None):
        """Help function for setting new Keywords in self._dictionary

        Args:
            keyword1 (str): keyword to set in self._dictionary
            keyword2 (str): keyword to set under keyword1
            keyword3 (str): keyword to set under keyword2
        """
        try:
            self._dictionary[keyword1]
        except KeyError:
            self._dictionary[keyword1] = {}
        if keyword2 is not None:
            try:
                self._dictionary[keyword1][keyword2]
            except KeyError:
                self._dictionary[keyword1][keyword2] = {}
        if keyword2 is not None and keyword3 is not None:
            try:
                self._dictionary[keyword1][keyword2][keyword3]
            except KeyError:
                self._dictionary[keyword1][keyword2][keyword3] = {}

    @staticmethod
    def read_max_angular_momentum(path):
        """Read maximum angular momentum from .skf file.
           See dftb.org/fileadmin/DFTB/public/misc/slakoformat.pdf
           for a detailed description of the Slater-Koster file format.

        Args:
            path (str): path to Slater-Koster file

        Returns:
            maxang (int/NoneType): max. angular momentum or
                None if extraction not possible
        """

        with open(path, 'r') as skf:
            line = skf.readline()
            if line[0] == '@':
                # Skip additional line for extended format
                line = skf.readline()

        # Replace any commas that may appear
        # (inconsistency in .skf files)
        line = line.replace(',', ' ').split()

        if len(line) == 3:
            # max. angular momentum specified:
            # extraction possible
            maxang = int(line[2]) - 1
        else:
            # max. angular momentum not specified
            # or wrong format: extraction not possible
            return None

        return maxang
