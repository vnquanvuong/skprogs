'''Defines supported xc-functionals.'''


import sktools.hsd as hsd
import sktools.hsd.converter as conv
import sktools.common as sc


class XCPBE0(sc.ClassDict):
    '''Globald PBE0 hybrid xc-functional.

    Attributes
    ----------
    alpha (float): fraction of the global exact HF exchange
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        alpha, child = query.getvalue(root, 'alpha', conv.float0,
                                      returnchild=True)
        if not 0.0 <= alpha <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid alpha CAM-parameter {:f}'.format(alpha),
                node=child)

        if not 0.0 <= alpha <= 1.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid global HFX portion alpha={:f}!\n'
                .format(alpha) +
                'Should satisfy 0.0 <= alpha <= 1.0', node=child)

        myself = cls()
        myself.type = 'pbe0'
        # dummy omega
        myself.omega = 1.0
        myself.alpha = alpha
        # dummy beta
        myself.beta = 0.0
        return myself


class XCB3LYP(sc.ClassDict):
    '''Globald B3LYP hybrid xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''

        myself = cls()
        myself.type = 'b3lyp'
        # dummy omega
        myself.omega = 1.0
        myself.alpha = 0.2
        myself.beta = 0.0
        return myself


class XCCAMYB3LYP(sc.ClassDict):
    '''Range-separated CAMY-B3LYP xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    alpha (float): fraction of the global exact HF exchange
    beta (float): determines (alpha + beta) fraction of long-range HF exchange
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        omega, child = query.getvalue(root, 'omega', conv.float0,
                                      returnchild=True)
        if omega <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid rs-parameter {:f}'.format(omega),
                node=child)

        alpha, child = query.getvalue(root, 'alpha', conv.float0,
                                      returnchild=True)

        beta, child = query.getvalue(root, 'beta', conv.float0,
                                     returnchild=True)

        if not alpha + beta > 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid CAM-parameter combination alpha={:f}, beta={:f}!\n'
                .format(alpha, beta) +
                'Should satisfy alpha + beta > 0.0', node=child)

        myself = cls()
        myself.type = 'camy-b3lyp'
        myself.omega = omega
        myself.alpha = alpha
        myself.beta = beta
        return myself


class XCCAMYPBEH(sc.ClassDict):
    '''Range-separated CAMY-PBEh xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    alpha (float): fraction of the global exact HF exchange
    beta (float): determines (alpha + beta) fraction of long-range HF exchange
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        omega, child = query.getvalue(root, 'omega', conv.float0,
                                      returnchild=True)
        if omega <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid rs-parameter {:f}'.format(omega),
                node=child)

        alpha, child = query.getvalue(root, 'alpha', conv.float0,
                                      returnchild=True)

        beta, child = query.getvalue(root, 'beta', conv.float0,
                                     returnchild=True)

        if not alpha + beta > 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid CAM-parameter combination alpha={:f}, beta={:f}!\n'
                .format(alpha, beta) +
                'Should satisfy alpha + beta > 0.0', node=child)

        myself = cls()
        myself.type = 'camy-pbeh'
        myself.omega = omega
        myself.alpha = alpha
        myself.beta = beta
        return myself


class XCLCYBNL(sc.ClassDict):
    '''Long-range corrected BNL xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        omega, child = query.getvalue(root, 'omega', conv.float0,
                                      returnchild=True)
        if omega <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid rs-parameter {:f}'.format(omega),
                node=child)

        myself = cls()
        myself.type = 'lcy-bnl'
        myself.omega = omega
        myself.alpha = 0.0
        myself.beta = 1.0
        return myself


class XCLCYPBE(sc.ClassDict):
    '''Long-range corrected PBE xc-functional.

    Attributes
    ----------
    omega (float): range-separation parameter
    '''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        omega, child = query.getvalue(root, 'omega', conv.float0,
                                      returnchild=True)
        if omega <= 0.0:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid rs-parameter {:f}'.format(omega),
                node=child)

        myself = cls()
        myself.type = 'lcy-pbe'
        myself.omega = omega
        myself.alpha = 0.0
        myself.beta = 1.0
        return myself


class XCLocal(sc.ClassDict):
    '''Local xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'local'
        return myself


class XCPBE(sc.ClassDict):
    '''Semi-local PBE xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'pbe'
        return myself


class XCBLYP(sc.ClassDict):
    '''Semi-local BLYP xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'blyp'
        return myself


class XCLDA(sc.ClassDict):
    '''Local LDA xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'lda'
        return myself


class XCr2SCAN(sc.ClassDict):
    '''Meta-GGA r2SCAN xc-functional.'''

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = 'r2scan'
        return myself


class _XCSimple(sc.ClassDict):
    '''Semilocal xc-functional with no parameters (type set by subclass _TYPE).'''
    _TYPE = None

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = cls._TYPE
        return myself


class _XCGlobalHybrid(sc.ClassDict):
    '''Global hybrid xc-functional with a fixed global HFX fraction (camAlpha), hardcoded in the
    Fortran codes and mirrored here in _ALPHA. It carries (omega, alpha, beta) = (1.0, camAlpha,
    0.0) so the SK-file RangeSep tag is written like the other hybrids; omega is a placeholder
    (there is no range separation, beta = 0).'''
    _TYPE = None
    _ALPHA = None

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myself = cls()
        myself.type = cls._TYPE
        myself.omega = 1.0
        myself.alpha = cls._ALPHA
        myself.beta = 0.0
        return myself


class _XCRangeSepErf(sc.ClassDict):
    '''erf range-separated hybrid (wB97X/wB97M) realised through a sum of Yukawa kernels (only the
    Yukawa kernel e^{-a r}/r has closed-form Slater-type-orbital two-electron integrals). The exact
    exchange is alpha*K_full + beta*K_erfLR with K_erfLR = sum_i c_i K_LRYukawa(beta_i*omega), the
    omega-independent {c_i, beta_i} hardcoded in the Fortran codes. M (myukawa, default 20) is the
    number of Yukawa terms approximating erf/erfc; VV10 nonlocal correlation is dropped (it is not
    representable in the SK-file framework -> pair with D3 for the -D3 form).

    Attributes
    ----------
    omega (float): range-separation parameter (0.30 for the wB97)
    alpha (float): short-range exact HF exchange fraction (camAlpha)
    beta (float): camBeta, so the long-range HF exchange fraction is alpha + beta
    myukawa (int): number of Yukawa terms M (1..30)
    '''
    _TYPE = None
    _OMEGA = 0.30
    _ALPHA = None
    _BETA = None

    @classmethod
    def fromhsd(cls, root, query):
        '''Creates instance from a HSD-node and with given query object.'''
        myukawa, child = query.getvalue(root, 'myukawa', conv.int0, defvalue=20,
                                        returnchild=True)
        if not 1 <= myukawa <= 30:
            raise hsd.HSDInvalidTagValueException(
                msg='Invalid number of Yukawa terms M={:d} (must be 1 <= M <= 30)'.format(myukawa),
                node=child)

        myself = cls()
        myself.type = cls._TYPE
        myself.omega = cls._OMEGA
        myself.alpha = cls._ALPHA
        myself.beta = cls._BETA
        myself.myukawa = myukawa
        return myself


class XCB97D(_XCSimple):    _TYPE = 'b97-d'
class XCM06L(_XCSimple):    _TYPE = 'm06-l'
class XCB97M(_XCSimple):    _TYPE = 'b97m'
class XCB972(_XCGlobalHybrid):   _TYPE = 'b97-2';   _ALPHA = 0.21
class XCB973(_XCGlobalHybrid):   _TYPE = 'b97-3';   _ALPHA = 0.269288
class XCr2SCANh(_XCGlobalHybrid): _TYPE = 'r2scanh'; _ALPHA = 0.10
class XCr2SCAN0(_XCGlobalHybrid): _TYPE = 'r2scan0'; _ALPHA = 0.25
class XCPW6B95(_XCGlobalHybrid):  _TYPE = 'pw6b95';  _ALPHA = 0.28
class XCMN15(_XCGlobalHybrid):    _TYPE = 'mn15';    _ALPHA = 0.44
class XCM062X(_XCGlobalHybrid):   _TYPE = 'm06-2x';  _ALPHA = 0.54
class XCwB97XV(_XCRangeSepErf):   _TYPE = 'wb97x-v'; _ALPHA = 0.167; _BETA = 0.833
class XCwB97MV(_XCRangeSepErf):   _TYPE = 'wb97m-v'; _ALPHA = 0.15;  _BETA = 0.85
# Phase-2 erf range-separated GGAs (omega/alpha/beta are hard-coded in the Fortran; the only HSD
# parameter is myukawa). Screened (HSE): alpha=-beta; LC: alpha=0,beta=1; CAM: general.
class XCHSE06(_XCRangeSepErf):    _TYPE = 'hse06';     _OMEGA = 0.11;  _ALPHA = 0.25;  _BETA = -0.25
class XCHSE12(_XCRangeSepErf):    _TYPE = 'hse12';     _OMEGA = 0.0978977840165; _ALPHA = 0.313; _BETA = -0.313
class XCLCwPBE(_XCRangeSepErf):   _TYPE = 'lc-wpbe';   _OMEGA = 0.40;  _ALPHA = 0.0;   _BETA = 1.0
class XCLCPBE(_XCRangeSepErf):    _TYPE = 'lc-pbe';    _OMEGA = 0.40;  _ALPHA = 0.0;   _BETA = 1.0
class XCLCBNL(_XCRangeSepErf):    _TYPE = 'lc-bnl';    _OMEGA = 0.33;  _ALPHA = 0.0;   _BETA = 1.0
class XCCAMB3LYP(_XCRangeSepErf): _TYPE = 'cam-b3lyp'; _OMEGA = 0.33;  _ALPHA = 0.19;  _BETA = 0.46
class XCCAMPBEh(_XCRangeSepErf):  _TYPE = 'cam-pbeh';  _OMEGA = 0.70;  _ALPHA = 1.0;   _BETA = -0.80
class XCWHPBE0(_XCRangeSepErf):   _TYPE = 'whpbe0';    _OMEGA = 0.20;  _ALPHA = 0.25;  _BETA = 0.25
# Phase-1 (non range-separated) functionals
class XCrevPBE(_XCSimple):        _TYPE = 'revpbe'
class XCRPBE(_XCSimple):          _TYPE = 'rpbe'
class XCTPSS(_XCSimple):          _TYPE = 'tpss'
class XCTASK(_XCSimple):          _TYPE = 'task'
class XCMN15L(_XCSimple):         _TYPE = 'mn15-l'
class XCrevPBE0(_XCGlobalHybrid):  _TYPE = 'revpbe0';  _ALPHA = 0.25
class XCTPSSh(_XCGlobalHybrid):    _TYPE = 'tpssh';    _ALPHA = 0.10
class XCr2SCAN50(_XCGlobalHybrid): _TYPE = 'r2scan50'; _ALPHA = 0.50
class XCM06(_XCGlobalHybrid):      _TYPE = 'm06';      _ALPHA = 0.27
class XCCF22D(_XCGlobalHybrid):    _TYPE = 'cf22d';    _ALPHA = 0.462806
# B97. Pure semilocal GGA:
class XCB973c(_XCSimple):         _TYPE = 'b97-3c'
# Global hybrid GGAs (camAlpha hard-coded in the Fortran):
class XCB97(_XCGlobalHybrid):     _TYPE = 'b97';   _ALPHA = 0.1943
class XCB971(_XCGlobalHybrid):    _TYPE = 'b97-1'; _ALPHA = 0.21
class XCB97K(_XCGlobalHybrid):    _TYPE = 'b97-k'; _ALPHA = 0.42
# Erf range-separated hybrid GGAs (omega/alpha/beta hard-coded in the Fortran; HSD param: myukawa):
class XCwB97(_XCRangeSepErf):     _TYPE = 'wb97';     _OMEGA = 0.40; _ALPHA = 0.0;      _BETA = 1.0
class XCwB97X(_XCRangeSepErf):    _TYPE = 'wb97x';    _OMEGA = 0.30; _ALPHA = 0.157706; _BETA = 0.842294
class XCwB97XD(_XCRangeSepErf):   _TYPE = 'wb97x-d';  _OMEGA = 0.20; _ALPHA = 0.222036; _BETA = 0.777964
class XCwB97XD3(_XCRangeSepErf):  _TYPE = 'wb97x-d3'; _OMEGA = 0.25; _ALPHA = 0.195728; _BETA = 0.804272
# O3LYP (global hybrid GGA, 11.61% HF; OPTX+LYP) and OPBE (pure GGA, OPTX + PBE):
class XCO3LYP(_XCGlobalHybrid):   _TYPE = 'o3lyp'; _ALPHA = 0.1161
class XCOPBE(_XCSimple):          _TYPE = 'opbe'


# Registered xc-functionals with corresponding HSD name as key:
XCFUNCTIONALS = {
    'lcy-bnl': XCLCYBNL,
    'lcy-pbe': XCLCYPBE,
    'local': XCLocal,
    'pbe': XCPBE,
    'lda': XCLDA,
    'blyp': XCBLYP,
    'pbe0': XCPBE0,
    'b3lyp': XCB3LYP,
    'camy-b3lyp': XCCAMYB3LYP,
    'camy-pbeh': XCCAMYPBEH,
    'r2scan': XCr2SCAN,
    'b97-d': XCB97D,
    'm06-l': XCM06L,
    'b97m': XCB97M,
    'b97-2': XCB972,
    'b97-3': XCB973,
    'r2scanh': XCr2SCANh,
    'r2scan0': XCr2SCAN0,
    'pw6b95': XCPW6B95,
    'mn15': XCMN15,
    'm06-2x': XCM062X,
    'wb97x-v': XCwB97XV,
    'wb97m-v': XCwB97MV,
    # wB97M-D3 shares wB97M-V's exchange parameters (the -D3 dispersion is a runtime DFTB+ add-on,
    # not encoded in the .skf), so it deliberately aliases the wB97M-V class.
    'wb97m-d3': XCwB97MV,
    'revpbe': XCrevPBE,
    'rpbe': XCRPBE,
    'tpss': XCTPSS,
    'task': XCTASK,
    'mn15-l': XCMN15L,
    'revpbe0': XCrevPBE0,
    'tpssh': XCTPSSh,
    'r2scan50': XCr2SCAN50,
    'm06': XCM06,
    'cf22d': XCCF22D,
    'hse06': XCHSE06,
    'hse12': XCHSE12,
    'lc-wpbe': XCLCwPBE,
    'lc-pbe': XCLCPBE,
    'lc-bnl': XCLCBNL,
    'cam-b3lyp': XCCAMB3LYP,
    'cam-pbeh': XCCAMPBEh,
    'whpbe0': XCWHPBE0,
    'b97-3c': XCB973c,
    'b97': XCB97,
    'b97-1': XCB971,
    'b97-k': XCB97K,
    'wb97': XCwB97,
    'wb97x': XCwB97X,
    'wb97x-d': XCwB97XD,
    'wb97x-d3': XCwB97XD3,
    'o3lyp': XCO3LYP,
    'opbe': XCOPBE
}
