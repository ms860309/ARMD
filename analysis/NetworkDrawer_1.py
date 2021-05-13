
from connect import db
import numpy as np
import os
import os.path as path
import sys
################################################################################
try:
    import cairocffi as cairo
except ImportError:
    try:
        import cairo
    except ImportError:
        cairo = None

sys.path.append(path.join(path.dirname(
    path.dirname(path.abspath(__file__))), 'database'))


def create_new_surface(file_format, target=None, width=1024, height=768):
    """
    Create a new surface of the specified `file_format`:
        "png" for :class:`ImageSurface`
        "svg" for :class:`SVGSurface`
        "pdf" for :class:`PDFSurface`
        "ps" for :class:`PSSurface`
    The surface will be written to the `target` parameter , which can be a
    path to save the surface to, or file-like object with a `write()` method.
    You can also optionally specify the `width` and `height` of the generated
    surface if you know what it is; otherwise a default size of 1024 by 768 is
    used.
    """
    file_format = file_format.lower()
    if file_format == 'png':
        surface = cairo.ImageSurface(
            cairo.FORMAT_ARGB32, int(width), int(height))
    elif file_format == 'svg':
        surface = cairo.SVGSurface(target, width, height)
    elif file_format == 'pdf':
        surface = cairo.PDFSurface(target, width, height)
    elif file_format == 'ps':
        surface = cairo.PSSurface(target, width, height)
    else:
        raise ValueError(
            'Invalid value "{0}" for type parameter; valid values are "png", "svg", "pdf", and "ps".'.format(type))
    return surface


class NetworkDrawer(object):
    """
    This class provides functionality for drawing the potential energy surface
    for a pressure-dependent reaction network using the Cairo 2D graphics
    engine. The most common use case is simply::

        NetworkDrawer().draw(network, file_format='png', path='network.png')

    where ``network`` is the :class:`Network` object to draw. You can also
    pass a dict of options to the constructor to affect how the network is
    drawn.
    """

    def __init__(self, options=None):
        self.options = {
            'structures': True,
            'fontFamily': 'sans',
            'fontSizeNormal': 3,
            'Eunits': 'kJ/mol',
            'padding': 16,
            'wellWidth': 64,
            'wellSpacing': 64,
            'Eslope': 1.5,
            'TSwidth': 16,
            'E0offset': 0.0,
        }
        if options:
            self.options.update(options)

        self.network = None
        self.dict_1 = {}
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None
        self.energy = {}

    def clear(self):
        self.network = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None

    def getEnergy(self):
        """
        PES energy is SCF energy, which could be calculate from qchem.
        """
        qm_collection = db['qm_calculate_center']
        for key, value in self.dict_1.items():
            query = {'path': value}
            targets = list(qm_collection.find(query))
            reactant_inchi_key = targets[0]['reactant_inchi_key']
            self.energy[reactant_inchi_key] = float(
                targets[0]['reactant_scf_energy'])*627.5095
        self.energy['AZHSSKPUVBVXLK-UHFFFAOYSA-N'] = -230.0389150086*627.5095
        self.energy['VHWYCFISAQVCCP-UHFFFAOYSA-N'] = -230.0385645234*627.5095
        self.energy['YRWZLJIAWSKYAM-UHFFFAOYSA-N'] = -230.0754621699*627.5095
        self.energy['YRUANFUOBUFKGQ-UHFFFAOYSA-N'] = -230.0231698005*627.5095
        self.energy['GJEMYANDWMDFCA-UHFFFAOYSA-N'] = -230.0689465192*627.5095
        self.energy['RZMLJGRDPXRMDB-UHFFFAOYSA-N'] = -230.0759236498*627.5095
        return self.energy

    def _get_energy_range(self):
        """
        Return the minimum and maximum energy on the potential energy
        surface.
        """

        energy = self.getEnergy()
        e0_min = energy[min(energy, key=energy.get)]
        e0_max = energy[max(energy, key=energy.get)]

        return e0_min, e0_max

    def _get_text_size(self, text, padding=2, file_format='pdf'):
        """

        """
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo

        # Use dummy surface to determine text extents
        surface = create_new_surface(file_format)
        cr = cairo.Context(surface)
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(text)
        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding

        return [0, 0, width, height]

    def _get_label_size(self, configuration, file_format='pdf'):
        """

        """
        bounding_rects = self._get_text_size(
            configuration, file_format=file_format)

        return bounding_rects

    def _get_total_dir(self, network):
        dirs = []
        for i in network:
            _dirs = self.network[i]
            for _dir in _dirs:
                if _dir not in dirs:
                    dirs.append(_dir)
        return dirs

    def _draw_text(self, text, cr, x0, y0, padding=2):
        """

        """
        cr.save()
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(text)
        cr.move_to(x0 - extents[0] - padding, y0 - extents[1] + padding)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(text)
        cr.restore()

        width = extents[2] + 2 * padding
        height = extents[3] + 2 * padding

        return [0, 0, width, height]

    def relative_height(self):
        """
        Use energy to get relative height and set reactant as 0
        """
        e_reactant = self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N']
        relative = {}
        for energy in self.energy:
            relative[energy] = self.energy[energy] - e_reactant
        return relative

    def draw(self, network, dict_2, file_format, path=None):
        # network
        self.dict_1 = dict_2
        self.network = network
        # draw parameter
        padding = self.options['padding']
        e_slope = self.options['Eslope']
        # determine the requiring height
        e_height = self._get_text_size('0.0', file_format=file_format)[3] + 6
        e0_min, e0_max = self._get_energy_range()
        height = (e0_max - e0_min) + 2 * padding + e_height
        length = {}
        for i in network:
            _length = len(network[i])
            length[i] = _length
        max_length_key = max(length, key=lambda k: length[k])
        max_length_value = length[max_length_key]
        width = 55 * max_length_value + 2 * padding - 40
        # Draw to the final surface
        surface = create_new_surface(
            file_format=file_format, target=path, width=width, height=height)
        cr = cairo.Context(surface)
        # draw
        # draw reactant horizontal line
        relative_height = self.relative_height()
        max_E_key = max(self.energy, key=lambda k: self.energy[k])
        max_E_value = self.energy[max_E_key]
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(1.0)
        cr.move_to(padding, max_E_value -
                   self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'] + 20)
        cr.line_to(padding + 10, max_E_value -
                   self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'] + 20)
        cr.save()
        # draw energy number and Energy unit
        """
        cr.move_to(padding, max_E_value - self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'] + 18)
        cr.set_font_size(3)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(str(round(self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'] - self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'], 3)))
        cr.restore()
        cr.save()
        """
        # unit
        cr.move_to(width - 30, height - 5)
        cr.set_font_size(3)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(str('Unit: kcal/mol'))
        cr.restore()
        cr.save()
        # draw dir label
        """
        cr.move_to(padding, max_E_value - self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'] + 23)
        cr.set_font_size(3)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text('CABDEMAGSHRORS-UHFFFAOYSA-N')
        cr.restore()
        cr.stroke()
        """
        #intermediate and product
        for reaction in network:
            for num, intermediate in enumerate(network[reaction][1:]):
                # draw intermediate and product horizontal line
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.set_line_width(0.5)
                cr.move_to(padding + 55 * (num+1), max_E_value -
                           self.energy[intermediate] + 20)
                cr.line_to(padding + 55 * (num+1) + 10,
                           max_E_value - self.energy[intermediate] + 20)
                cr.save()
                # draw energy number

                cr.move_to(padding + 55 * (num+2), max_E_value -
                           self.energy[intermediate] + 18)
                cr.set_font_size(3)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(str(round(
                    self.energy[intermediate] - self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'], 3)))
                cr.restore()
                cr.save()
                """
                #draw dir label
                cr.move_to(padding + 55 * (num+1), max_E_value - self.energy[network[reaction][num+1]] + 23)
                cr.set_font_size(1)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(intermediate + '('+ str(round(self.energy[intermediate] - self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'], 3)) +')')
                cr.restore()
                cr.save()
                """
                if len(network[reaction]) == 2:
                    cr.move_to(padding + 10, max_E_value -
                               self.energy['CABDEMAGSHRORS-UHFFFAOYSA-N'] + 20)
                    cr.line_to(padding + 55 * (num+1), max_E_value -
                               self.energy[intermediate] + 20)
                    cr.restore()
                    cr.save()
                elif self.network[reaction][-1] != intermediate:
                    cr.move_to(padding + 55 * (num+1) + 10,
                               max_E_value - self.energy[intermediate] + 20)
                    cr.line_to(padding + 55 * (num+2), max_E_value -
                               self.energy[network[reaction][num+2]] + 20)
                    cr.restore()
                    cr.save()
                cr.stroke()


reaction_collection = db['reactions']
query = {'$and':
         [{"unique":
           {"$in":
            ['new one']
            }
           }, {'for_debug':
               {'$nin':
                ['from same']}
               }]}

targets = list(reaction_collection.find(query))

dict_1 = {}
dict_2 = {}
for idx, i in enumerate(targets):
    _path = i['path']
    rxn = i['reaction']
    gen = i['generations']
    name = 'reaction_{}'.format(idx)
    dict_1[name] = _path
    dict_2[name] = [rxn, gen]

dict_2 = {
    'reaction_0': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'RZKDJWZLOVYVCV-UHFFFAOYSA-N', 'AZHSSKPUVBVXLK-UHFFFAOYSA-N'],
    'reaction_1': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'MUOXUGRTQNZNJE-UHFFFAOYSA-N', 'AZHSSKPUVBVXLK-UHFFFAOYSA-N'],
    'reaction_2': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'LYCAIKOWRPUZTN-UHFFFAOYSA-N', 'UPUZLJDFULGKCC-UHFFFAOYSA-N'],
    'reaction_3': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'LYCAIKOWRPUZTN-UHFFFAOYSA-N', 'VHWYCFISAQVCCP-UHFFFAOYSA-N'],
    'reaction_4': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'UPUZLJDFULGKCC-UHFFFAOYSA-N', 'YRWZLJIAWSKYAM-UHFFFAOYSA-N'],
    'reaction_5': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'UPUZLJDFULGKCC-UHFFFAOYSA-N', 'AZHSSKPUVBVXLK-UHFFFAOYSA-N'],
    'reaction_6': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'UPUZLJDFULGKCC-UHFFFAOYSA-N', 'YRUANFUOBUFKGQ-UHFFFAOYSA-N'],
    'reaction_7': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'LYCAIKOWRPUZTN-UHFFFAOYSA-N', 'REHUGJYJIZPQAV-UHFFFAOYSA-N', 'GJEMYANDWMDFCA-UHFFFAOYSA-N'],
    'reaction_8': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'LYCAIKOWRPUZTN-UHFFFAOYSA-N', 'REHUGJYJIZPQAV-UHFFFAOYSA-N', 'RZMLJGRDPXRMDB-UHFFFAOYSA-N'],
    'reaction_9': ['CABDEMAGSHRORS-UHFFFAOYSA-N', 'LYCAIKOWRPUZTN-UHFFFAOYSA-N', 'REHUGJYJIZPQAV-UHFFFAOYSA-N', 'VHWYCFISAQVCCP-UHFFFAOYSA-N']
}
NetworkDrawer().draw(dict_2, dict_1, 'pdf', os.path.join(os.getcwd(), 'network.pdf'))
