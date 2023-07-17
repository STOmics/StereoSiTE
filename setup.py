from  setuptools import setup, find_packages

setup(
    name='StereoSiTE',
    version='1.0',
    author='LiuXing',
    author_email='liuxing2@genomics.cn',
    description=('Analysis spatial transcriptomics data'),
    license='MIT License',
    keywords='spatial cell interaction intensity',

    packages=['stereosite'], #需要打包的目录列表

    include_package_data=True,
    platforms='any',
    #需要安装的依赖包
    install_requires = [
        'anndata>=0.8.0',
        'scanpy>=1.9.1',
        'squidpy>=1.1.2',
    ],

    entry_points = {
        'console_scripts': [
            'scii = stereosite.SCII:main'
            'deconcolution = stereosite.CN.deconvolution:main'
        ]
    },

    #程序的所属分类列表
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Bioinfomatics'
        'License :: OSI Approved :: MIT License'
    ],
)