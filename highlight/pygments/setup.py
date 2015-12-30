from setuptools import setup
 
setup(
    name='Playmol Pygments Lexer',
    version='0.0.1',
    description=__doc__,
    author='Charlles Abreu',
    packages=['playmol_lexer'],
    entry_points='''[pygments.lexers]
playmollexer = playmol_lexer:PlaymolLexer
'''
)

