from  setuptools import setup, find_packages

with open("PYPI.rst", "r") as f:
    long_description = f.read()

setup(
    name='stereosite',
    version='2.2.0',
    author='LiuXing',
    author_email='liuxing2@genomics.cn',
    description=('Analysis spatial transcriptomics data'),
    long_description=long_description,
    license='GPL-3 License',
    keywords='spatial cell interaction intensity',
    url="https://github.com/STOmics/StereoSiTE",

    packages=find_packages(), #['stereosite'], #需要打包的目录列表

    include_package_data=True,
    platforms='any',
    #需要安装的依赖包
    install_requires = [
        'anndata>=0.8.0',
        'scanpy>=1.9.1',
        'squidpy>=1.1.2',
        'decoupler>=1.4.0',
        'pydeseq2>=0.3.6',
        'networkx>=3.1',
        'tensorly>=0.8.1',
        'scikit-learn>=1.2.1',
        'torch>=1.11.0',
        'igraph>=0.10.4',
        'pycirclize>=1.1.0',
        'cell2location==0.1.3'
    ],
)
