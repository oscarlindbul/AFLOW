import json
import subprocess
import os

class XtalFinder:

    def __init__(self, aflow_executable='aflow'):
        self.aflow_executable = aflow_executable

    def aflow_command(self, cmd, input=None):
        try:
            return subprocess.check_output(
                self.aflow_executable + cmd,
                shell=True,
                universal_newlines=True,
                input=input
            )
        except subprocess.CalledProcessError:
            raise OSError('aflow executable not found: ' + self.aflow_executable)


    def check_path(self, path):
        if type(path) == str:
            if os.path.exists(path):
                return os.path.realpath(path).replace(" ","\\ ")
            else:
                raise OSError(filename + ' not found')
        elif type(path) == list:
            for p in path:
                if not os.path.exists(p):
                    raise OSError(p + ' not found')
            return ','.join(list(map(os.path.realpath, path))).replace(" ","\\ ")
        else:
            raise TypeError('The path to input file/files/directory must be a string or a list.')


    def get_prototype_label(self, input_file, options=None):
        fpath = self.check_path(input_file)
        command = ' --prototype'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def compare_materials(self, input_files, options=None):
        fpath = self.check_path(input_files)
        command = ' --compare_materials=' + fpath
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_materials_directory(self, directory, options=None):
        fpath = self.check_path(directory)
        command = ' --compare_materials -D ' + fpath
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_materials_file(self, filename, options=None):
        fpath = self.check_path(filename)
        command = ' --compare_materials -F=' + fpath
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_structures(self, input_files, options=None):
        fpath = self.check_path(input_files)
        command = ' --compare_structures=' + fpath
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json
    
    def compare_structures_string(self, *strings, options=None):
        command = ' --compare_structures'
        output = ''

        total_string = ""
        for s in strings:
            total_string += "[VASP_POSCAR_MODE_EXPLICIT]START\n"
            total_string += s
            total_string += "[VASP_POSCAR_MODE_EXPLICIT]STOP\n"

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet',
            input=total_string
        )

        res_json = json.loads(output)
        return res_json

    def compare_structures_directory(self, directory, options=None):
        fpath = self.check_path(directory)
        command = ' --compare_structures -D ' + fpath
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_structures_file(self, filename, options=None):
        fpath = self.check_path(filename)
        command = ' --compare_structures -F=' + fpath
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )
        print(output)

        res_json = json.loads(output)
        return res_json

    def compare2database(self, input_file, options=None):
        fpath = self.check_path(input_file)
        command = ' --compare2database'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def compare2prototypes(self, input_file, options=None):
        fpath = self.check_path(input_file)
        command = ' --compare2prototypes'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def get_isopointal_prototypes(self, input_file, options=None):
        fpath = self.check_path(input_file)
        command = ' --isopointal_prototype'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def get_unique_atom_decorations(self, input_file, options=None):
        fpath = self.check_path(input_file)
        command = ' --unique_atom_decorations'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json < ' + fpath
        )

        res_json = json.loads(output)
        return res_json
