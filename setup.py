from setuptools import setup, find_packages
from mockinbird import __version__

setup(
    name='mockinbird',
    version=__version__,
    description='PAR-CLIP data analysis pipeline',
    author='Christian Roth, Phillipp Torkler',
    author_email='christan.roth@mpibpc.mpg.de',
    license='GPLv3',
    entry_points={
        'console_scripts': [

            # general
            'mockinbird = mockinbird.scripts.mockinbird:main',
            'mockinbird-preprocess = mockinbird.scripts.preprocess:main',
            'mockinbird-postprocess = mockinbird.scripts.postprocess:main',

            # preprocessing
            'mb-remove-duplicates = mockinbird.scripts.removePCRduplicates:run',
            'mb-adapter-clipper = mockinbird.scripts.clipper53:main',
            'mb-bam-postprocess = mockinbird.scripts.bam_postprocessing:main',
            'mb-softclip-analyzer = mockinbird.scripts.clipped_seq:main',
            'mb-bsfinder = mockinbird.scripts.bsfinder:run',
            'mb-normalize = mockinbird.scripts.normalize:run',
            'mb-cap-occupancy = mockinbird.scripts.convert2quantile:run',
            'mb-table2fasta = mockinbird.utils.table2fasta:main',
            'mb-upgrade-table = mockinbird.utils.update_table:main',

            # postprocessing
            'mb-plot-metagene-nobs = mockinbird.plots.makeCenterBothEnds:run',
            'mb-plot-metagene = mockinbird.plots.makeCenterBothEnds_bs:main',
            'mb-plot-kmer-enrichment = mockinbird.plots.makeKmerPerPosition:run',
            'mb-generate-negative-set = mockinbird.scripts.makeNegSets:run',
            'mb-plot-kmer-logodds = mockinbird.plots.makeKmerLogOdds:run',
            'mb-xxmotif = mockinbird.scripts.xxmotif:run',
            'mb-plot-transition-frequencies = mockinbird.plots.makeNucleotideProbabilities:run',
            'mb-plot-heatmap = mockinbird.plots.makeHeatMap:run',
            'mb-plot-heatmap-small = mockinbird.plots.makeHeatMapSmall:run',
            'mb-filter-sites = mockinbird.scripts.filter_sites:main',
            'mb-annotate-table = mockinbird.scripts.utils.annotate_table:main',
        ]
    },
    packages=find_packages(),
    include_package_data=True,
    test_suite='nose.collector',
    zip_safe=False
)
