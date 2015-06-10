from setuptools import setup, find_packages

setup(
    name="selectionTools",
    version="1.1",
    packages=['selection_pipeline','selection_pipeline.tests'],
    test_suite='selection_pipeline.tests.test_selection_pipeline',
    author="James Boocock",
    author_email="smilefreak@gmx.com",
    description="Selection Pipeline for VCF Data",
    license="MIT",
    keywords="iHS ehh selection evolution",
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'ancestral_annotation = selection_pipeline.aa_annotate:main',
            'selection_pipeline = selection_pipeline.selection_pipeline:main',
            'multipop_selection_pipeline = selection_pipeline.multipipeline:main',
            'haps_to_hapmap = selection_pipeline.haps_to_hapmap:main',
            'haps_filters = selection_pipeline.haps_filters:main',
            'haps_interpolate = selection_pipeline.haps_interpolate:main'
        ]
    },
    url="github.com/smilefreak/MerrimanSelectionPipeline",
    use_2to3=True,
    include_package_data=True,
    package_data = {
        '' : ['*.haps','*.cfg','*.vcf','*.sample','*.ped',
              '*.map',"*.fa","*.ids"]
    }

)
