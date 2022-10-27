from setuptools import setup

setup(name='kmcalc',
      version='1.0.0',
      author=['Sergio Garcia','David Dooley'],
      packages=['kmcalc'],
      entry_points={
          'console_scripts': [
              'kmcalc = kmcalc.__main__:cli',
              'kmcalcgui = kmcalc.gui:cli',
          ]
      },
      )