from distutils.core import setup

setup (name = "Reconstructing_Macrocomplex",
        version = "1.0",
        description = "This program models macrocomplex structures of biomolecules formed by proteins. The macrocomplex is build from pairwaise interactions",
        author = "Irene Agustí Barea, Marta Espinosa Camarena",
        author_email = 'ireneagustibarea@gmail.com, espinosacamarenamarta@gmail.com',
        packages = ['Reconstructing_Macrocomplex'],
        requires = ['Bio', 'sys', 'os', 'argparse','logging','random','string','gzip','shutil'],
        scripts = ['Reconstructing_Macrocomplex/main_reconstructing_macrocomplex.py', 'Reconstructing_Macrocomplex/functions_reconstructing_macrocomplex.py']
)
