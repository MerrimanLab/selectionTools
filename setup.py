from setuptools import setup, find_packages

setup(
    name="selectionTools",
    version="0.9",
    packages=find_packages(),
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
            'multi_population = selection_pipeline.multipipeline:main'
        ]
    },
    url="github.com/smilefreak/MerrimanSelectionPipeline",
    use_2to3=True

)
