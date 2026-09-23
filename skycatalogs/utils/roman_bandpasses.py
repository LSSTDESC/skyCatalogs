"""Lazy access to optional Roman effective-area bandpasses."""

import galsim
import numpy as np


# Preserve the established SkyCatalog flux-column names.
ROMAN_BANDS = ('W146', 'R062', 'Z087', 'Y106', 'J129', 'H158', 'F184',
               'K213')
_TECHNICAL_NAMES = {
    'W146': 'F146',
    'R062': 'F062',
    'Z087': 'F087',
    'Y106': 'F106',
    'J129': 'F129',
    'H158': 'F158',
    'F184': 'F184',
    'K213': 'F213',
    'SNPrism': 'Prism',
    'Grism_1stOrder': 'Grism_1stOrder',
    'Grism_0thOrder': 'Grism_0thOrder',
}

_bandpasses = None
_version = None


def _load_roman_bandpasses(include_all_bands=False):
    """Return nominal Roman curves and their provenance string."""
    global _bandpasses, _version
    if _bandpasses is None:
        try:
            import roman_technical_information as roman_info
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                'Roman bandpasses require the optional dependency installed '
                'by `pip install skyCatalogs[roman]`') from exc

        effective_area_dir = (
            roman_info.PACKAGEDIR
            / 'WideFieldInstrument/Imaging/EffectiveAreas')
        resources = sorted(
            item for item in effective_area_dir.iterdir()
            if item.name.startswith('Roman_effarea_v8_SCA')
            and item.name.endswith('.ecsv'))
        if len(resources) != 18:
            raise RuntimeError(
                'Expected 18 Roman imaging effective-area SCA tables; found '
                f'{len(resources)}')

        column_names = None
        wavelength_micron = None
        sca_values = []
        for resource in resources:
            with resource.open('r') as stream:
                for line in stream:
                    if line.lstrip().startswith('Wave,'):
                        names = tuple(item.strip() for item in line.split(','))
                        break
                else:
                    raise RuntimeError(
                        f'No data header found in {resource.name}')
                values = np.loadtxt(stream, delimiter=',')
            if column_names is None:
                column_names = names
                wavelength_micron = values[:, 0]
            elif names != column_names or not np.array_equal(
                    values[:, 0], wavelength_micron):
                raise RuntimeError(
                    'Roman effective-area SCA tables do not share one grid')
            sca_values.append(values[:, 1:])

        mean_effective_area = np.mean(sca_values, axis=0)
        column_index = {name: index - 1
                        for index, name in enumerate(column_names) if index}
        wavelength_nm = wavelength_micron * 1000.0
        _bandpasses = {}
        for name, technical_name in _TECHNICAL_NAMES.items():
            throughput = mean_effective_area[:, column_index[technical_name]]
            peak = np.max(throughput)
            if peak <= 0:
                raise RuntimeError(
                    f'Roman effective-area curve {technical_name} is empty')
            table = galsim.LookupTable(
                wavelength_nm, throughput / peak, interpolant='linear')
            bandpass = galsim.Bandpass(table, wave_type='nm')
            _bandpasses[name] = bandpass.truncate(
                relative_throughput=1.e-4).thin().withZeropoint('AB')
        _version = (
            f'roman-technical-information-{roman_info.__version__};'
            'mean-SCA01-18;relative-throughput-cut=1e-4')

    names = _TECHNICAL_NAMES if include_all_bands else ROMAN_BANDS
    return {name: _bandpasses[name] for name in names}, _version


def load_roman_bandpasses(include_all_bands=False):
    """Return Roman bandpasses."""
    return _load_roman_bandpasses(include_all_bands)[0]


def load_roman_bandpasses_version():
    """Return provenance for the Roman bandpasses."""
    return _load_roman_bandpasses()[1]
