import os
import shelve
import subprocess as subproc
import numpy as np
from sktools import hsd
import sktools.hsd.converter as conv
import sktools.common as sc
from sktools import twocenter_grids
from sktools import radial_grid


# Category-grouped, contiguous IDs (must match sktwocnt/lib/xcfunctionals.f90; = slateratom - 1).
SUPPORTED_FUNCTIONALS = {
    # LDA
    'lda' : 1,
    # GGA (pure)
    'pbe' : 2, 'revpbe' : 3, 'rpbe' : 4, 'blyp' : 5, 'b97-d' : 6, 'b97-3c' : 7, 'opbe' : 8,
    # meta-GGA (pure)
    'r2scan' : 9, 'b97m' : 10, 'tpss' : 11, 'task' : 12, 'm06-l' : 13, 'mn15-l' : 14,
    # global hybrid + GGA
    'pbe0' : 15, 'b3lyp' : 16, 'b97' : 17, 'b97-1' : 18, 'b97-2' : 19, 'b97-3' : 20,
    'b97-k' : 21, 'revpbe0' : 22, 'o3lyp' : 23,
    # global hybrid + meta-GGA
    'r2scanh' : 24, 'r2scan0' : 25, 'r2scan50' : 26, 'pw6b95' : 27, 'tpssh' : 28,
    'm06' : 29, 'm06-2x' : 30, 'mn15' : 31, 'cf22d' : 32,
    # range-separated + GGA (pure LC)
    'lcy-pbe' : 33, 'lcy-bnl' : 34, 'lc-wpbe' : 35, 'lc-pbe' : 36, 'lc-bnl' : 37, 'wb97' : 38,
    # range-separated + global hybrid + GGA
    'hse06' : 39, 'hse12' : 40, 'camy-b3lyp' : 41, 'camy-pbeh' : 42, 'cam-b3lyp' : 43,
    'cam-pbeh' : 44, 'whpbe0' : 45, 'wb97x' : 46, 'wb97x-d' : 47, 'wb97x-d3' : 48, 'wb97x-v' : 49,
    # range-separated + global hybrid + meta-GGA
    'wb97m-v' : 50}

INPUT_FILE = "sktwocnt.in"
STDOUT_FILE = "output"
BASISFUNCTION_FILE = "basisfuncs.dbm"
DEFAULT_BINARY = "sktwocnt"


class SktwocntSettings(sc.ClassDict):
    """Specific settings for sktwocnt program.

    Attributes
    ----------
    integrationpoints : int, int
        Two integers representing the nr. of points for radial and angular
        integration.
    """

    def __init__(self, integrationpoints):
        super().__init__()
        self.integrationpoints = integrationpoints

    @classmethod
    def fromhsd(cls, node, query):
        """Generate the object from HSD tree"""
        integrationpoints, child = query.getvalue(
            node, "integrationpoints", conv.int1, returnchild=True)
        if len(integrationpoints) != 2:
            raise hsd.HSDInvalidTagValueException(
                "Two integration point parameters must be specified", child)
        return cls(integrationpoints)

    def __eq__(self, other):
        if not isinstance(other, SktwocntSettings):
            return False
        if self.integrationpoints != other.integrationpoints:
            return False
        return True


class Sktwocnt:

    def __init__(self, workdir):
        self._workdir = workdir

    def set_input(self, settings, superpos, functional, grid, atom1data,
                  atom2data=None):
        myinput = SktwocntInput(settings, superpos, functional, grid, atom1data,
                                atom2data)
        myinput.write(self._workdir)

    def run(self, binary=DEFAULT_BINARY):
        runner = SktwocntCalculation(binary, self._workdir)
        runner.run()

    def get_result(self):
        result = SktwocntResult(self._workdir)
        return result


class SktwocntInput:

    _INTERACTION_FROM_NTYPES = {
        1: "homo",
        2: "hetero",
    }

    _POTENTIAL_SUPERPOS = "potential"

    def __init__(self, settings, superpos, functional, grid, atom1data,
                 atom2data=None):
        self._settings = settings
        self._atom1data = atom1data
        self._hetero = atom2data is not None
        if self._hetero:
            self._atom2data = atom2data
        else:
            self._atom2data = self._atom1data
        self._check_superposition(superpos)
        self._densitysuperpos = (superpos == sc.SUPERPOSITION_DENSITY)
        self._check_functional(functional.type)
        self._functional = functional
        self._check_grid(grid)
        self._grid = grid

    @staticmethod
    def _check_superposition(superpos):
        if superpos not in \
        [sc.SUPERPOSITION_POTENTIAL, sc.SUPERPOSITION_DENSITY]:
            msg = "Sktwocnt: Invalid superposition type"
            sc.SkgenException(msg)

    @staticmethod
    def _check_functional(functional):
        if functional not in SUPPORTED_FUNCTIONALS:
            raise sc.SkgenException("Invalid functional type")

    @staticmethod
    def _check_grid(grid):
        if not isinstance(grid, twocenter_grids.EquidistantGrid):
            msg = "Sktwocnt can only handle equidistant grids"
            raise sc.SkgenException(msg)

    def write(self, workdir):
        atomfiles1 = self._store_atomdata(workdir, self._atom1data, 1)
        if self._hetero:
            atomfiles2 = self._store_atomdata(workdir, self._atom2data, 2)
        else:
            atomfiles2 = None
        self._store_twocnt_input(workdir, atomfiles1, atomfiles2)
        self._store_basisfunctions(workdir)

    def _store_atomdata(self, workdir, atomdata, iatom):
        atomfiles = sc.ClassDict()
        atomfiles.wavefuncs = self._store_wavefuncs(workdir, atomdata.wavefuncs,
                                                    iatom)
        atomfiles.potential = self._store_potentials(workdir,
                                                     atomdata.potentials, iatom)
        atomfiles.density = self._store_density(workdir, atomdata.density,
                                                iatom)
        xcn = self._functional.type
        if xcn in ('lcy-bnl', 'lcy-pbe', 'pbe0', 'b3lyp', 'camy-b3lyp',
                   'camy-pbeh', 'b97-2', 'b97-3', 'r2scanh', 'r2scan0',
                   'pw6b95', 'mn15', 'm06-2x', 'wb97x-v', 'wb97m-v',
                   'revpbe0', 'tpssh', 'r2scan50', 'm06', 'cf22d',
                   'hse06', 'hse12', 'lc-wpbe', 'lc-pbe', 'lc-bnl',
                   'cam-b3lyp', 'cam-pbeh', 'whpbe0',
                   'b97', 'b97-1', 'b97-k', 'o3lyp',
                   'wb97', 'wb97x', 'wb97x-d', 'wb97x-d3'):
            atomfiles.dens_wavefuncs = self._store_dens_wavefuncs(
                workdir, atomdata.dens_wavefuncs, iatom)
        atomfiles.occshells = atomdata.occshells
        return atomfiles

    @staticmethod
    def _store_wavefuncs(workdir, wavefuncs, iatom):
        wavefuncfiles = []
        for nn, ll, wfc012 in wavefuncs:
            fname = "wave{:d}_{:d}{:s}.dat".format(iatom, nn,
                                                   sc.ANGMOM_TO_SHELL[ll])
            wfc012.tofile(os.path.join(workdir, fname))
            wavefuncfiles.append((nn, ll, fname))
        return wavefuncfiles

    @staticmethod
    def _store_dens_wavefuncs(workdir, wavefuncs, iatom):
        wavefuncfiles = []
        for nn, ll, wfc012 in wavefuncs:
            fname = "dens_wave{:d}_{:d}{:s}.dat".format(iatom, nn,
                                                        sc.ANGMOM_TO_SHELL[ll])
            wfc012.tofile(os.path.join(workdir, fname))
            wavefuncfiles.append((nn, ll, fname))
        return wavefuncfiles

    @staticmethod
    def _store_potentials(workdir, potentials, iatom):
        fname = "potentials{:d}.dat".format(iatom)
        # Vxc up and down should be equivalent, twocnt reads only one.
        newdata = potentials.data.take((radial_grid.VNUC, radial_grid.VHARTREE,
                                        radial_grid.VXCUP), axis=1)
        newgriddata = radial_grid.GridData(potentials.grid, newdata)
        newgriddata.tofile(os.path.join(workdir, fname))
        return fname

    @staticmethod
    def _store_density(workdir, density, iatom):
        fname = "density{:d}.dat".format(iatom)
        density.tofile(os.path.join(workdir, fname))
        return fname

    def _store_basisfunctions(self, workdir):
        config = shelve.open(
            os.path.join(workdir, BASISFUNCTION_FILE), "n")
        config["basis1"] = [(nn, ll) for nn, ll, wfc012
                            in self._atom1data.wavefuncs]
        config["basis2"] = [(nn, ll) for nn, ll, wfc012
                            in self._atom2data.wavefuncs]
        config.close()

    def _store_twocnt_input(self, workdir, atomfiles1, atomfiles2=None):
        fp = open(os.path.join(workdir, INPUT_FILE), "w")
        self._write_twocnt_header(fp)
        self._write_twocnt_gridinfo(fp)
        self._write_twocnt_integration_parameters(fp)
        self._write_twocnt_atom_block(fp, atomfiles1)
        if self._hetero:
            self._write_twocnt_atom_block(fp, atomfiles2)
        fp.close()

    def _write_twocnt_header(self, fp):
        if self._densitysuperpos:
            superposname = 'density'
        else:
            superposname = 'potential'

        xcfkey = self._functional.type
        ixc = SUPPORTED_FUNCTIONALS[xcfkey]
        fp.write('{} {} {}\n'.format('hetero' if self._hetero else 'homo',
                                     superposname, ixc))

    def _write_twocnt_gridinfo(self, fp):
        '''Writes integration grid info.'''

        # long-range corrected functionals
        if self._functional.type in ('lcy-bnl', 'lcy-pbe'):
            # hardcoded parameters for the Becke integration,
            # -> should probably be moved to skdef.hsd
            becke = '2000 194 11 1.0'
            fp.write("{:f}\n".format(self._functional.omega))
            fp.write("{:s}\n".format(becke))
        # B3LYP and other global hybrids with fixed HFX fraction (alpha hard-coded in sktwocnt):
        # only the Becke grid is written
        elif self._functional.type in ('b3lyp', 'b97-2', 'b97-3', 'r2scanh',
                                        'r2scan0', 'pw6b95', 'mn15', 'm06-2x',
                                        'revpbe0', 'tpssh', 'r2scan50', 'm06', 'cf22d',
                                        'b97', 'b97-1', 'b97-k', 'o3lyp'):
            # hardcoded parameters for the Becke integration,
            # -> should probably be moved to skdef.hsd
            becke = '2000 194 11 1.0'
            fp.write("{:s}\n".format(becke))
        # PBE0
        elif self._functional.type == 'pbe0':
            becke = '2000 194 11 1.0'
            fp.write("{:f}\n".format(self._functional.alpha))
            fp.write("{:s}\n".format(becke))
        # CAM functionals
        elif self._functional.type in ('camy-b3lyp', 'camy-pbeh'):
            becke = '2000 194 11 1.0'
            fp.write("{:f} {:f} {:f}\n".format(self._functional.omega,
                                               self._functional.alpha,
                                               self._functional.beta))
            fp.write("{:s}\n".format(becke))
        # erf range-separated hybrids (wB97X-V/wB97M-V, Phase-2 screened/LC/CAM GGAs, B97 wB97):
        # number of Yukawa terms M, then the Becke grid (omega, camAlpha, camBeta hard-coded in sktwocnt)
        elif self._functional.type in ('wb97x-v', 'wb97m-v', 'hse06', 'hse12', 'lc-wpbe', 'lc-pbe',
                                        'lc-bnl', 'cam-b3lyp', 'cam-pbeh', 'whpbe0',
                                        'wb97', 'wb97x', 'wb97x-d', 'wb97x-d3'):
            becke = '2000 194 11 1.0'
            fp.write("{:d}\n".format(self._functional.myukawa))
            fp.write("{:s}\n".format(becke))

        fp.write("{:f} {:f} {:e} {:f}\n".format(
            self._grid.gridstart, self._grid.gridseparation,
            self._grid.tolerance, self._grid.maxdistance))

    def _write_twocnt_integration_parameters(self, fp):
        fp.write("{:d} {:d}\n".format(*self._settings.integrationpoints))

    def _write_twocnt_atom_block(self, fp, atomfiles):
        if self._functional.type in ('lcy-bnl', 'lcy-pbe', 'pbe0', 'b3lyp',
                                     'camy-b3lyp', 'camy-pbeh', 'b97-2', 'b97-3',
                                     'r2scanh', 'r2scan0', 'pw6b95', 'mn15', 'm06-2x',
                                     'wb97x-v', 'wb97m-v', 'revpbe0', 'tpssh', 'r2scan50',
                                     'm06', 'cf22d', 'hse06', 'hse12', 'lc-wpbe', 'lc-pbe',
                                     'lc-bnl', 'cam-b3lyp', 'cam-pbeh', 'whpbe0',
                                     'b97', 'b97-1', 'b97-k', 'o3lyp',
                                     'wb97', 'wb97x', 'wb97x-d', 'wb97x-d3'):
            fp.write("{:d} {:d}\n".format(len(atomfiles.wavefuncs),
                                          len(atomfiles.dens_wavefuncs)))
        else:
            fp.write("{:d}\n".format(len(atomfiles.wavefuncs)))

        for nn, ll, wavefuncfile in atomfiles.wavefuncs:
            fp.write("'{}' {:d}\n".format(wavefuncfile, ll))

        if self._functional.type in ('lcy-bnl', 'lcy-pbe', 'pbe0', 'b3lyp',
                                     'camy-b3lyp', 'camy-pbeh', 'b97-2', 'b97-3',
                                     'r2scanh', 'r2scan0', 'pw6b95', 'mn15', 'm06-2x',
                                     'wb97x-v', 'wb97m-v', 'revpbe0', 'tpssh', 'r2scan50',
                                     'm06', 'cf22d', 'hse06', 'hse12', 'lc-wpbe', 'lc-pbe',
                                     'lc-bnl', 'cam-b3lyp', 'cam-pbeh', 'whpbe0',
                                     'b97', 'b97-1', 'b97-k', 'o3lyp',
                                     'wb97', 'wb97x', 'wb97x-d', 'wb97x-d3'):
            occdict = {}
            for xx in atomfiles.occshells:
                occdict[xx[0]] = xx[1]

            for nn, ll, dens_wavefuncfile in atomfiles.dens_wavefuncs:
                fp.write("'{}' {:d} {:f}\n"
                         .format(dens_wavefuncfile, ll, occdict[(nn, ll)]))

        fp.write("'{}'\n".format(atomfiles.potential))
        if self._densitysuperpos:
            fp.write("'{}'\n".format(atomfiles.density))
        else:
            fp.write("'{}'\n".format("nostart"))


class SktwocntCalculation:

    def __init__(self, binary, workdir):
        self._binary = binary
        self._workdir = workdir

    def run(self):
        fpin = open(os.path.join(self._workdir, INPUT_FILE), "r")
        fpout = open(os.path.join(self._workdir, STDOUT_FILE), "w")
        proc = subproc.Popen(self._binary.split(), cwd=self._workdir,
                             stdin=fpin, stdout=fpout, stderr=subproc.STDOUT)
        proc.wait()
        fpin.close()
        fpout.close()


class SktwocntResult:

    def __init__(self, workdir):
        basis1, basis2 = self._read_basis(workdir)
        self._integmap = self._create_integral_mapping(basis1, basis2)
        ninteg = len(self._integmap)
        self._skham = self._read_sktable(
            os.path.join(workdir, "at1-at2.ham.dat"), ninteg)
        self._skover = self._read_sktable(
            os.path.join(workdir, "at1-at2.over.dat"), ninteg)

    @staticmethod
    def _read_basis(workdir):
        config = shelve.open(os.path.join(workdir, BASISFUNCTION_FILE), "r")
        basis1 = list(config["basis1"])
        basis2 = list(config["basis2"])
        config.close()
        return basis1, basis2

    @staticmethod
    def _create_integral_mapping(basis1, basis2):
        ninteg = 0
        integmap = {}
        for n1, l1 in basis1:
            for n2, l2 in basis2:
                for mm in range(min(l1, l2) + 1):
                    ninteg += 1
                    integmap[(n1, l1, n2, l2, mm)] = ninteg
        return integmap

    @staticmethod
    def _read_sktable(fname, ninteg):
        fp = open(fname, "r")
        nline = int(fp.readline())
        # noinspection PyNoneFunctionAssignment,PyTypeChecker
        tmp = np.fromfile(fp, dtype=float, count=ninteg * nline, sep=" ")
        tmp.shape = (nline, ninteg)
        return tmp

    def get_hamiltonian(self):
        return self._skham

    def get_overlap(self):
        return self._skover
