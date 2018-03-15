import os
import time
from functools import partial
from mockinbird.utils import pipeline as pl
from mockinbird.utils import config_validation as cv


class BamStatisticsModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        r_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('gff_file', cv.Annot(str, default='', converter=r_conv)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline

        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']

        prefix = general_cfg['prefix']

        bam_file = pipeline.get_curfile(fmt='bam')

        stat_file = os.path.join(output_dir, prefix + '_stat.json')
        cmd = [
            'mb-create-bam-statistics',
            '%r' % bam_file,
            '%r' % stat_file,
        ]

        if cfg['gff_file']:
            cmd.append('--gff_file %r' % cfg['gff_file'])

        self._cmds.append(cmd)
        self._intermed_files.append(stat_file)
        pipeline.upd_curfile(fmt='stat_file', filepath=stat_file)



class Bam2MockinbirdModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        rw_conv = partial(cv.rel_file_rw_validator, cfg_path=pipeline.cfg_path)
        r_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('mock_bam', cv.Annot(converter=r_conv)),
            ('genome_fasta', cv.Annot(converter=r_conv)),
            ('bed_regions', cv.Annot(converter=r_conv, default='', warn_if_missing=False)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline
        general_cfg = pipeline.get_config('general')
        output_dir = general_cfg['output_dir']
        read_cfg = pipeline.get_config('reads')
        prefix = general_cfg['prefix']

        factor_bam = pipeline.get_curfile(fmt='bam')
        mock_bam = cfg['mock_bam']

        site_freqs = os.path.join(output_dir, prefix + '.site_freqs')
        binding_sites = os.path.join(output_dir, prefix + '.sites')

        cmd = [
            'bam2mockinbird',
            '%r' % factor_bam,
            '%r' % mock_bam,
            '%r' % cfg['genome_fasta'],
            '%r' % binding_sites,
            '%r' % site_freqs,
            '--transition_from %s' % read_cfg['reference_nucleotide'],
            '--transition_to %s' % read_cfg['mutation_nucleotide'],
        ]
        if cfg['bed_regions']:
            cmd.append('-l %r' % cfg['bed_regions'])

        self._cmds.append(cmd)

        pipeline.upd_curfile(fmt='sites', filepath=binding_sites)
        pipeline.upd_curfile(fmt='site_freqs', filepath=site_freqs)


class LearnMockModule(pl.CmdPipelineModule):

    def __init__(self, pipeline):
        rw_conv = partial(cv.rel_file_rw_validator, cfg_path=pipeline.cfg_path)
        r_conv = partial(cv.rel_file_r_validator, cfg_path=pipeline.cfg_path)
        cfg_fmt = [
            ('mock_model', cv.Annot(converter=rw_conv)),
            ('mock_statistics', cv.Annot(converter=r_conv)),
            ('n_mixture_components', cv.Annot(converter=cv.nonneg_integer)),
            ('em_iterations', cv.Annot(default=250, converter=cv.nonneg_integer)),
        ]
        super().__init__(pipeline, cfg_req=cfg_fmt)

    def prepare(self, cfg):
        super().prepare(cfg)
        pipeline = self._pipeline

        full_table = pipeline.get_curfile(fmt='site_freqs')
        mock_model = cfg['mock_model']

        mock_parent = os.path.dirname(mock_model)

        learn_cmd = [
            'mb-learn-mock',
            '%r' % full_table,
            cfg['n_mixture_components'],
            cfg['mock_statistics'],
            mock_parent,
            '--n_iterations %s' % cfg['em_iterations'],

        ]
        if not os.path.exists(mock_model):
            self._cmds.append(learn_cmd)

        pipeline.upd_curfile(fmt='mock_model', filepath=mock_model)
