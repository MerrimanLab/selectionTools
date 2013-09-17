from setuptools import setup, find_packages

setup(
    name="MerrimanSelectionPipeline",
    version="0.3",
    packages=find_packages(),
            
    author= "James Boocock",
    author_email="smilefreak@gmx.com",
    description="Selection Pipeline for VCF Data",
    license = "MIT",
    keywords="iHS ehh selection evolution",
    zip_safe=False,
    entry_points ={
        'console_scripts': [
            'ancestral_annotation = selection_pipeline.aa_annotate:main',
            'selection_pipeline = selection_pipeline.fullprocess:main',
        ]
    },
    url="github.com/smilefreak/MerrimanSelectionPipeline"


)
