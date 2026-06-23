############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, Jul 2013           #
############################################################


# glob 是 Python 标准库，用来按照通配符查找文件，例如 '*.png' 表示所有 png 图片。
import glob
# os 是 Python 标准库，用来和操作系统交互，例如拼接路径、创建目录、切换目录、读取环境变量。
import os
# shutil 是 Python 标准库，用来做高级文件操作，例如复制文件 copy2()、移动文件 move()。
import shutil
# time 是 Python 标准库，用来记录时间；这里用于统计某个处理步骤花了多久。
import time

# numpy 是常用的科学计算库；这里用 np.sum() 等函数对数组像元进行统计。
import numpy as np

# mintpy 是当前项目的主包；mintpy.__file__ 可以告诉我们 MintPy 安装目录在哪里。
import mintpy
# RAMP_LIST 是允许的相位坡度类型列表；cluster 管理并行计算；sensor 处理卫星传感器名称。
from mintpy.objects import RAMP_LIST, cluster, sensor
# readfile/writefile 负责读写 HDF5/模板等文件；utils as ut 是把 utils 模块起短名 ut，方便后面调用。
from mintpy.utils import readfile, utils as ut, writefile

# 说明：本文件负责把 MintPy 小基线时序分析的各个命令行模块串成完整工作流。
# “工作流”可以理解为：按固定顺序调用 load_data、modify_network、ifgram_inversion 等工具。
# comment out the workflow import in favor of individual imports for speed and robustness
# import mintpy.workflow  # dynamic import of modules for smallbaselineApp


##########################################################################
def get_the_latest_default_template_file(work_dir):
    """Get the latest version of default template file.
    If an obsolete file exists in the working directory, the existing option values are kept.
    """
    # 获取软件自带的最新默认模板，以及当前工作目录中的默认模板路径。
    # get the latest and current versions of the default template files
    # os.path.dirname(mintpy.__file__) 得到 mintpy 包所在目录。
    # os.path.join(...) 用系统认可的路径分隔符拼出 defaults/smallbaselineApp.cfg 的完整路径。
    lfile = os.path.join(os.path.dirname(mintpy.__file__), 'defaults/smallbaselineApp.cfg')
    # work_dir 是用户要运行项目的目录；cfile 是这个目录里的默认配置文件路径。
    cfile = os.path.join(work_dir, 'smallbaselineApp.cfg')

    # os.path.isfile(cfile) 判断 cfile 这个路径是否存在且确实是一个文件。
    if not os.path.isfile(cfile):
        # 工作目录中没有模板时，直接复制一份默认模板。
        print(f'copy default template file {lfile} to work directory')
        # shutil.copy2(src, dst) 复制文件，同时尽量保留原文件的时间戳等元信息。
        shutil.copy2(lfile, work_dir)
    else:
        # 工作目录已有模板时，检查是否缺少新版模板中的配置项。
        #cfile is obsolete if any key is missing
        # readfile.read_template() 是 MintPy 的函数，会把 cfg 模板读成 Python 字典 dict。
        ldict = readfile.read_template(lfile)
        cdict = readfile.read_template(cfile)
        # any([...]) 表示只要列表里有一个 True，整体就是 True；这里用来判断旧模板是否缺 key。
        if any([key not in cdict.keys() for key in ldict.keys()]):
            print('obsolete default template detected, update to the latest version.')
            shutil.copy2(lfile, work_dir)

            # 更新模板结构后，保留旧模板中用户已经设置过的参数值。
            #keep the existing option value from obsolete template file
            # ut.update_template_file() 会把 cdict 中已有的用户配置写回新版模板文件。
            ut.update_template_file(cfile, cdict)

    # 返回当前工作目录里最终可用的模板文件路径。
    return cfile


##########################################################################
class TimeSeriesAnalysis:
    """ Routine processing workflow for time series analysis of small baseline InSAR stacks
    """

    def __init__(self, customTemplateFile=None, workDir=None):
        # 记录用户模板、工作目录以及启动程序时所在目录，便于流程结束后返回。
        # __init__ 是 Python 类的初始化方法；创建 TimeSeriesAnalysis(...) 对象时会自动运行。
        # self 表示“当前这个对象本身”；self.xxx 是把变量保存到对象里，供其它方法继续使用。
        self.customTemplateFile = customTemplateFile
        # os.path.abspath() 把相对路径转换为绝对路径，避免后续切换目录后找不到文件。
        self.workDir = os.path.abspath(workDir)
        # os.getcwd() 表示 get current working directory，即获取当前终端所在目录。
        self.cwd = os.path.abspath(os.getcwd())

    def open(self):
        """The starting point of the workflow. It runs every time.
        It 1) grab project name if given
           2) go to work directory
           3) get and read template(s) options
        """

        #1. Get projectName
        # projectName 默认没有；如果用户提供了模板，就用模板文件名当项目名。
        self.projectName = None
        if self.customTemplateFile:
            # os.path.basename() 取文件名，os.path.splitext() 把文件名和扩展名分开。
            self.projectName = os.path.splitext(os.path.basename(self.customTemplateFile))[0]
            print('Project name:', self.projectName)

        #2. Go to work directory
        # os.makedirs(..., exist_ok=True) 创建目录；目录已存在时不会报错。
        os.makedirs(self.workDir, exist_ok=True)
        # os.chdir() 切换 Python 当前工作目录；后面相对路径都会以 self.workDir 为起点。
        os.chdir(self.workDir)
        print('Go to work directory:', self.workDir)

        #3. Read templates
        # 调用本文件上面定义的函数，确保工作目录里有最新版 smallbaselineApp.cfg。
        self.templateFile = get_the_latest_default_template_file(self.workDir)
        # 调用当前类自己的 _read_template() 方法，读取并合并模板参数。
        self._read_template()


    def _read_template(self):
        # 读取并合并用户模板与默认模板，是后续所有处理步骤的参数来源。
        # 方法名前面的下划线 _ 表示这个方法主要给类内部使用，不是给用户直接调用的入口。
        # 1) update default template
        self.customTemplate = None
        if self.customTemplateFile:
            # customTemplateFile --> customTemplate
            print('read custom template file:', self.customTemplateFile)
            # 把用户模板文件读成字典，例如 {'mintpy.geocode': True, ...}。
            cdict = readfile.read_template(self.customTemplateFile)

            # 将常见的简写/布尔写法规范化，减少模板输入格式差异造成的问题。
            # correct some loose type errors
            standardValues = {
                'def':'auto', 'default':'auto',
                'y':'yes', 'on':'yes', 'true':'yes',
                'n':'no', 'off':'no', 'false':'no',
            }
            for key, value in cdict.items():
                # cdict.items() 会逐个返回字典中的 key 和 value。
                if value in standardValues.keys():
                    cdict[key] = standardValues[value]

            for key in ['mintpy.deramp', 'mintpy.troposphericDelay.method']:
                if key in cdict.keys():
                    # lower() 转小写，replace('-', '_') 把连字符替换成下划线，统一方法名格式。
                    cdict[key] = cdict[key].lower().replace('-', '_')

            # 兼容旧模板里直接使用 processor 的写法。
            if 'processor' in cdict.keys():
                cdict['mintpy.load.processor'] = cdict['processor']

            # 以下偏移元数据只在 load_data.py 中临时使用，合并模板后不再需要保留。
            # these metadata are used in load_data.py only, not needed afterwards
            # (in order to manually add extra offset when the lookup table is shifted)
            # (seen in ROI_PAC product sometimes)
            for key in ['SUBSET_XMIN', 'SUBSET_YMIN']:
                if key in cdict.keys():
                    # pop() 会从字典中删除这个 key，并返回被删除的值。
                    cdict.pop(key)

            # dict(cdict) 复制一份字典，保存为对象属性，后面写元数据时还会用到。
            self.customTemplate = dict(cdict)

            # 用用户模板中的配置覆盖默认模板对应项，生成本次运行使用的模板文件。
            # customTemplate --> templateFile
            print('update default template based on input custom template')
            self.templateFile = ut.update_template_file(self.templateFile, self.customTemplate)

        # 将模板备份到 inputs 和 pic 目录，便于追踪本次处理使用的参数。
        # 2) backup custom/default template file in inputs/pic folder
        flen = len(os.path.basename(self.templateFile))
        if self.customTemplateFile:
            # max(a, b) 取较大值；这里用于对齐打印信息中的文件名宽度。
            flen = max(flen, len(os.path.basename(self.customTemplateFile)))

        for backup_dirname in ['inputs', 'pic']:
            # 循环创建 inputs 和 pic 两个备份目录。
            backup_dir = os.path.join(self.workDir, backup_dirname)
            # create directory
            os.makedirs(backup_dir, exist_ok=True)

            # back up to the directory
            # 列表推导式 [x for x in ... if x] 会过滤掉 None 或空字符串。
            tfiles = [x for x in [self.customTemplateFile, self.templateFile] if x]
            for tfile in tfiles:
                # out_file 是备份后的目标文件路径。
                out_file = os.path.join(backup_dir, os.path.basename(tfile))
                # ut.run_or_skip() 是 MintPy 的工具函数，用来判断目标文件是否需要重新生成/复制。
                if ut.run_or_skip(out_file, in_file=tfile, readable=False, print_msg=False) == 'run':
                    shutil.copy2(tfile, backup_dir)
                    print('copy {f:<{l}} to {d:<8} directory for backup.'.format(
                        f=os.path.basename(tfile), l=flen, d=os.path.basename(backup_dir)))

        # 3) read default template file
        print('read default template file:', self.templateFile)
        # self.template 保存最终生效的完整模板配置，后面每个处理步骤都会读取它。
        self.template = readfile.read_template(self.templateFile)
        # check_template_auto_value() 会把模板里的 auto/yes/no 等文本进一步转换成程序需要的值。
        self.template = ut.check_template_auto_value(self.template)

        # 若用户要求输出 HDF-EOS5/KMZ，则需要先开启地理编码。
        # correct some loose setup conflicts
        if self.template['mintpy.geocode'] is False:
            for key in ['mintpy.save.hdfEos5', 'mintpy.save.kmz']:
                if self.template[key] is True:
                    self.template['mintpy.geocode'] = True
                    print(f'Turn ON mintpy.geocode in order to run {key}.')
                    break


    def run_load_data(self, step_name):
        """Load InSAR stacks into HDF5 files in ./inputs folder.
        It 1) copy auxiliary files into work directory (for Unvi of Miami only)
           2) load all interferograms stack files into mintpy/inputs directory.
           3) check loading result
           4) add custom metadata (optional, for HDF-EOS5 format only)
        """
        # load_data 是工作流入口步骤：把外部 InSAR 产品整理成 MintPy HDF5 输入。
        # 1) copy aux files (optional)
        # self._copy_aux_file() 调用当前类中的辅助方法，尝试复制一些可选的外部辅助文件。
        self._copy_aux_file()

        # 2) loading data
        # compose list of input arguments
        # instead of using command line then split
        # to support path with whitespace
        # iargs 是“参数列表”，等价于用户在终端里输入 load_data.py 后面的那些参数。
        # 用列表而不是一整串字符串，可以正确处理路径中包含空格的情况。
        iargs = ['--template', self.templateFile]
        if self.customTemplateFile:
            # 如果用户提供了自定义模板，也把它传给 load_data.py。
            iargs += [self.customTemplateFile]
        if self.projectName:
            # 如果能从模板名得到项目名，就通过 --project 参数传给 load_data.py。
            iargs += ['--project', self.projectName]

        # run command line
        print('\nload_data.py', ' '.join(iargs))
        # import mintpy.cli.load_data 表示导入 MintPy 里的 load_data 命令模块。
        # 这里不是开启新的终端进程，而是在 Python 代码里直接调用它的 main() 函数。
        import mintpy.cli.load_data
        # main(iargs) 会按 iargs 里的参数执行 load_data 的功能。
        mintpy.cli.load_data.main(iargs)

        # come back to working directory
        # load_data 执行过程中可能切换目录，所以这里强制切回本次项目工作目录。
        os.chdir(self.workDir)

        # 3) check loading result
        # ut.check_loaded_dataset() 会检查 inputs 目录中的关键 HDF5 文件是否已经正确生成。
        # [:4] 表示只取返回结果的前 4 个：干涉图栈、几何文件、查找表、电离层文件。
        stack_file, geom_file, _, ion_file = ut.check_loaded_dataset(self.workDir, print_msg=True)[:4]

        # 若提供了自定义模板，把其中的元数据补写到加载后的数据文件中。
        # 4) add custom metadata (optional)
        if self.customTemplateFile:
            # use ut.add_attribute() instead of add_attribute.py because of
            # better control of special metadata, such as SUBSET_X/YMIN
            # msg 是日志文字，f'...' 是 Python 的格式化字符串，可以把变量值嵌入文本。
            msg = f'updating metadata based on custom template file {os.path.basename(self.customTemplateFile)}'
            for fname in [stack_file, ion_file]:  #, geom_file]:
                if fname:
                    print(f'{msg} for file: {os.path.basename(fname)}')
                    # ut.add_attribute() 会把 self.customTemplate 中的键值对写入 HDF5 文件属性。
                    ut.add_attribute(fname, self.customTemplate)


    def _copy_aux_file(self):
        # 针对 University of Miami 环境的辅助文件拷贝逻辑，普通运行环境通常不会触发。
        # for Univ of Miami
        # os.getenv('SCRATCHDIR') 用来读取名为 SCRATCHDIR 的环境变量；没有设置时返回 None。
        if os.getenv('SCRATCHDIR') and self.projectName:
            # proj_dir 是外部临时/共享目录中当前项目的位置。
            proj_dir = os.path.join(os.getenv('SCRATCHDIR'), self.projectName)
            flist = ['PROCESS/unavco_attributes.txt',
                     'PROCESS/bl_list.txt',
                     'SLC/summary*slc.jpg']
            # ut.get_file_list() 会把通配符路径展开成真实文件列表；abspath=True 表示返回绝对路径。
            flist = ut.get_file_list([os.path.join(proj_dir, i) for i in flist], abspath=True)
            for fname in flist:
                if ut.run_or_skip(os.path.basename(fname), in_file=fname, readable=False) == 'run':
                    # 把辅助文件复制到当前工作目录。
                    shutil.copy2(fname, self.workDir)
                    print(f'copy {os.path.basename(fname)} to work directory')


    def run_network_modification(self, step_name):
        """Modify network of interferograms before the network inversion."""
        # 在反演前裁剪或调整干涉网络，并准备网络图所需的输入。
        # check the existence of ifgramStack.h5
        # stack_file 是干涉图栈文件；geom_file 是几何信息文件。
        stack_file, geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[:2]
        # coherenceSpatialAvg.txt 是空间平均相干性文本文件，plot_network 会用到。
        coh_txt = os.path.join(self.workDir, 'coherenceSpatialAvg.txt')
        # net_fig 同时检查工作目录和 pic 目录下是否已有 network.pdf。
        net_fig = [os.path.join(self.workDir, i, 'network.pdf') for i in ['', 'pic']]
        net_fig = [x for x in net_fig if os.path.isfile(x)]

        # 1) output waterMask.h5 to simplify the detection/use of waterMask
        water_mask_file = os.path.join(self.workDir, 'waterMask.h5')
        # readfile.get_dataset_list() 会列出 HDF5 文件里有哪些数据集，例如 waterMask、latitude。
        if 'waterMask' in readfile.get_dataset_list(geom_file):
            if ut.run_or_skip(out_file=water_mask_file, in_file=geom_file) == 'run':
                print(f'generate {water_mask_file} from {geom_file} for conveniency')
                # readfile.read() 从 HDF5 文件读取指定数据集；返回值通常是数据数组和属性字典。
                water_mask, atr = readfile.read(geom_file, datasetName='waterMask')

                # 经纬度为 0 的像元通常是几何文件中的无效区域，也同步标记为水体/无效。
                # ignore no-data pixels in geometry files
                ds_name_list = readfile.get_dataset_list(geom_file)
                for ds_name in ['latitude','longitude']:
                    if ds_name in ds_name_list:
                        print(f'set pixels with 0 in {ds_name} to 0 in waterMask')
                        ds = readfile.read(geom_file, datasetName=ds_name)[0]
                        # water_mask[ds == 0] = 0 是 numpy 数组筛选写法：把经纬度为 0 的位置设为 0。
                        water_mask[ds == 0] = 0

                atr['FILE_TYPE'] = 'waterMask'
                # writefile.write() 把数组和元数据写成新的 HDF5 文件。
                writefile.write(water_mask, out_file=water_mask_file, metadata=atr)

        # 2) modify network
        iargs = [stack_file, '-t', self.templateFile]
        print('\nmodify_network.py', ' '.join(iargs))
        # 调用 modify_network.py，根据模板配置删除低质量干涉图或修改网络连接。
        import mintpy.cli.modify_network
        mintpy.cli.modify_network.main(iargs)

        # 3) plot network
        iargs = [stack_file, '-t', self.templateFile, '--nodisplay']

        # 相位栈和偏移栈使用不同的数据集作为网络质量指标。
        dsNames = readfile.get_dataset_list(stack_file)
        # any(...) 判断数据集名称中是否包含 phase/offset，从而决定绘图时用 coherence 还是 offsetSNR。
        if any('phase' in i.lower() for i in dsNames):
            iargs += ['-d', 'coherence', '-v', '0.2', '1.0']
        elif any('offset' in i.lower() for i in dsNames):
            iargs += ['-d', 'offsetSNR', '-v', '0', '20']

        print('\nplot_network.py', ' '.join(iargs))

        # run
        if self.template['mintpy.plot']:
            if ut.run_or_skip(out_file=net_fig,
                              in_file=[stack_file, coh_txt, self.templateFile],
                              readable=False) == 'run':
                # 调用 plot_network.py 生成干涉网络图。
                import mintpy.cli.plot_network
                mintpy.cli.plot_network.main(iargs)
        else:
            print('mintpy.plot is turned OFF, skip plotting network.')


    def generate_ifgram_aux_file(self):
        """Generate auxiliary files from ifgramStack file"""
        # 生成参考点选择和后续反演常用的干涉图辅助文件。
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[0]
        dsNames = readfile.get_dataset_list(stack_file)
        # 下面三个文件分别保存连接分量掩膜、平均空间相干性、平均空间信噪比。
        mask_file = os.path.join(self.workDir, 'maskConnComp.h5')
        coh_file = os.path.join(self.workDir, 'avgSpatialCoh.h5')
        snr_file = os.path.join(self.workDir, 'avgSpatialSNR.h5')

        # 1) generate mask file from the common connected components
        if any('phase' in i.lower() for i in dsNames):
            iargs = [stack_file, '--nonzero', '-o', mask_file, '--update']
            print('\ngenerate_mask.py', ' '.join(iargs))
            # generate_mask.py 根据非零像元生成掩膜文件。
            import mintpy.cli.generate_mask
            mintpy.cli.generate_mask.main(iargs)

        # 2) generate average spatial coherence
        # temporal_average.py 虽然名字是时间平均，也可用于对干涉图栈中的 coherence/offsetSNR 求平均。
        if any('phase' in i.lower() for i in dsNames):
            iargs = [stack_file, '--dataset', 'coherence', '-o', coh_file, '--update']
        elif any('offset' in i.lower() for i in dsNames):
            iargs = [stack_file, '--dataset', 'offsetSNR', '-o', snr_file, '--update']
        print('\ntemporal_average.py', ' '.join(iargs))
        import mintpy.cli.temporal_average
        mintpy.cli.temporal_average.main(iargs)


    def run_reference_point(self, step_name):
        """Select reference point.
        It 1) generate mask file from common conn comp
           2) generate average spatial coherence and its mask
           3) add REF_X/Y and/or REF_LAT/LON attribute to stack file
        """
        # 选择参考像元，并将参考像元信息写入干涉图栈属性。
        # 1-2) aux files: maskConnComp and avgSpatialCoh
        # 先确保参考点选择所需的辅助文件已经生成。
        self.generate_ifgram_aux_file()

        # 3) add REF_X/Y(/LAT/LON) of the reference point
        # lookup_file 是雷达坐标到地理坐标的查找表；没有查找表时可能为 None。
        stack_file, _, lookup_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[:3]
        coh_file = os.path.join(self.workDir, 'avgSpatialCoh.h5')

        iargs = [stack_file, '-t', self.templateFile, '-c', coh_file]
        if lookup_file is not None:
            iargs += ['--lookup', lookup_file]
        print('\nreference_point.py', ' '.join(iargs))
        # reference_point.py 会根据模板和相干性文件选择参考点，并写入 REF_X/REF_Y 等属性。
        import mintpy.cli.reference_point
        mintpy.cli.reference_point.main(iargs)


    def run_quick_overview(self, step_name):
        """A quick overview on the interferogram stack for:
            1) avgPhaseVelocity.h5: possible ground deformation through interferogram stacking
            2) numTriNonzeroIntAmbiguity.h5: phase unwrapping errors
                through the integer ambiguity of phase closure
        """
        # 快速生成两个诊断产品，用于初步检查形变信号和解缠错误。
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[0]

        # 1) stack interferograms
        pha_vel_file = os.path.join(self.workDir, 'avgPhaseVelocity.h5')
        iargs = [stack_file, '--dataset', 'unwrapPhase', '-o', pha_vel_file, '--update']
        print('\ntemporal_average.py', ' '.join(iargs))
        # 对 unwrapPhase 做时间平均，得到一个快速查看整体形变趋势的文件。
        import mintpy.cli.temporal_average
        mintpy.cli.temporal_average.main(iargs)

        # 2) calculate the number of interferogram triplets with non-zero integer ambiguity
        water_mask_file = os.path.join(self.workDir, 'waterMask.h5')
        iargs = [stack_file, '--water-mask', water_mask_file, '--action', 'calculate', '--update']
        print('\nunwrap_error_phase_closure.py', ' '.join(iargs))
        # unwrap_error_phase_closure.py 这里用于计算相位闭合诊断量，而不是直接改正。
        import mintpy.cli.unwrap_error_phase_closure
        mintpy.cli.unwrap_error_phase_closure.main(iargs)


    def run_unwrap_error_correction(self, step_name):
        """Correct phase-unwrapping errors"""
        # 根据模板中指定的方法执行解缠误差改正；未指定方法则跳过。
        # self.template 是模板字典；用方括号读取某个配置项的值。
        method = self.template['mintpy.unwrapError.method']
        if not method:
            print('phase-unwrapping error correction is OFF.')
            return

        # check required input files
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[0]
        mask_file = os.path.join(self.workDir, 'maskConnComp.h5')

        iargs_bridge = [stack_file, '--template', self.templateFile, '--update']
        iargs_closure = iargs_bridge + ['--cc-mask', mask_file]

        # 先导入两个可能会用到的 MintPy 命令模块，后面根据 method 选择调用哪个 main()。
        import mintpy.cli.unwrap_error_bridging
        import mintpy.cli.unwrap_error_phase_closure
        if method == 'bridging':
            # 基于连接分量桥接的解缠误差修正。
            print('\nunwrap_error_bridging.py', ' '.join(iargs_bridge))
            mintpy.cli.unwrap_error_bridging.main(iargs_bridge)

        elif method == 'phase_closure':
            # 基于相位闭合的解缠误差修正。
            print('\nunwrap_error_phase_closure.py', ' '.join(iargs_closure))
            mintpy.cli.unwrap_error_phase_closure.main(iargs_closure)

        elif method == 'bridging+phase_closure':
            # 组合方法：先桥接修正，再进行相位闭合修正。
            # -i 指定输入数据集名，-o 指定输出数据集名，避免覆盖原始 unwrapPhase。
            iargs_bridge += ['-i', 'unwrapPhase',
                             '-o', 'unwrapPhase_bridging']
            print('\nunwrap_error_bridging.py', ' '.join(iargs_bridge))
            mintpy.cli.unwrap_error_bridging.main(iargs_bridge)

            iargs_closure += ['-i', 'unwrapPhase_bridging',
                              '-o', 'unwrapPhase_bridging_phaseClosure']
            print('\nunwrap_error_phase_closure.py', ' '.join(iargs_closure))
            mintpy.cli.unwrap_error_phase_closure.main(iargs_closure)

        else:
            # raise ValueError 表示输入值不合法，程序会停止并显示错误信息。
            raise ValueError(f'un-recognized method: {method}')


    def run_network_inversion(self, step_name):
        """Invert network of interferograms for raw phase time-series.
        1) network inversion --> timeseries.h5, temporalCoherence.h5, numInvIfgram.h5
        2) temporalCoherence.h5 --> maskTempCoh.h5
        """
        # 将干涉网络反演为原始时序，并生成可靠像元掩膜。
        # check the existence of ifgramStack.h5
        stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[0]

        # 1) invert ifgramStack for time-series
        iargs = [stack_file, '-t', self.templateFile, '--update']
        print('\nifgram_inversion.py', ' '.join(iargs))
        # ifgram_inversion.py 是核心反演步骤，把干涉图网络转换为位移时间序列。
        import mintpy.cli.ifgram_inversion
        mintpy.cli.ifgram_inversion.main(iargs)

        # 2) get reliable pixel mask: maskTempCoh.h5
        # 调用本类自己的方法，根据 temporalCoherence.h5 生成可靠像元掩膜。
        self.generate_temporal_coherence_mask()


    def generate_temporal_coherence_mask(self):
        """Generate reliable pixel mask from temporal coherence"""
        # 根据时间相干性阈值生成可靠像元掩膜，供后续时序分析使用。
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        # temporalCoherence.h5 是网络反演产生的时间相干性文件。
        tcoh_file = os.path.join(self.workDir, 'temporalCoherence.h5')
        mask_file = os.path.join(self.workDir, 'maskTempCoh.h5')
        tcoh_min = self.template['mintpy.networkInversion.minTempCoh']

        # compose list of arguments
        # -m 是最小阈值；-o 是输出文件名。
        iargs = [tcoh_file, '-m', tcoh_min, '-o', mask_file]
        # exclude pixels in shadow if shadowMask dataset is available
        if (self.template['mintpy.networkInversion.shadowMask'] is True
                and 'shadowMask' in readfile.get_dataset_list(geom_file)):
            # 如果几何文件中有 shadowMask，就把阴影区域排除在可靠像元之外。
            iargs += ['--base', geom_file, '--base-dataset', 'shadowMask', '--base-value', '1']
        print('\ngenerate_mask.py', ' '.join(iargs))

        # update 模式只在输出过期或关键配置变化时重新生成掩膜。
        # update mode: run only if:
        # 1) output file exists and newer than input file, AND
        # 2) all config keys are the same
        config_keys = [f'mintpy.networkInversion.{i}' for i in ['minTempCoh','shadowMask']]
        print('update mode: ON')
        flag = 'skip'
        if ut.run_or_skip(out_file=mask_file, in_file=tcoh_file, print_msg=False) == 'run':
            flag = 'run'
        else:
            print(f'1) output file: {mask_file} already exists and newer than input file: {tcoh_file}')
            # readfile.read_attribute() 只读取文件属性，不读取大数组数据，速度更快。
            atr = readfile.read_attribute(mask_file)
            if any(str(self.template[i]) != atr.get(i, 'False') for i in config_keys):
                flag = 'run'
                print(f'2) NOT all key configuration parameters are the same: {config_keys}')
            else:
                print(f'2) all key configuration parameters are the same: {config_keys}')
        print(f'run or skip: {flag}')

        if flag == 'run':
            import mintpy.cli.generate_mask
            mintpy.cli.generate_mask.main(iargs)

            # 将关键配置写入掩膜文件属性，便于下次判断是否需要重新运行。
            # update config_keys
            atr = {}
            for key in config_keys:
                atr[key] = self.template[key]
            ut.add_attribute(mask_file, atr)

        # check number of pixels selected in mask file for following analysis
        # readfile.read(mask_file)[0] 读取掩膜数组；!= 0 找出有效像元；np.sum() 统计数量。
        num_pixel = np.sum(readfile.read(mask_file)[0] != 0.)
        print(f'number of reliable pixels: {num_pixel}')

        min_num_pixel = float(self.template['mintpy.networkInversion.minNumPixel'])
        if num_pixel < min_num_pixel:
            # 可靠像元过少时，后续分析结果不可信，因此提前停止并生成诊断图。
            msg = f"Not enough reliable pixels (minimum of {int(min_num_pixel)}). "
            msg += "Try the following:\n"
            msg += "1) Check the reference pixel and make sure it's not in areas with unwrapping errors\n"
            msg += "2) Check the network and make sure it's fully connected without subsets"
            print(f'ERROR: {msg}')

            # populate the pic folder to facilate the trouble shooting
            # 如果可靠像元太少，先自动绘图，方便用户排查参考点、网络连接等问题。
            self.plot_result(print_aux=False)

            # terminate the program
            raise RuntimeError(msg)


    @staticmethod
    def get_timeseries_filename(template, work_dir='./'):
        """Get input/output time-series filename for each step
        Parameters: template : dict, content of smallbaselineApp.cfg
        Returns:    steps    : dict of dicts, input/output filenames for each step
        """
        # 根据模板中开启的改正项，推导每一步的输入/输出时序文件名。
        work_dir = os.path.abspath(work_dir)
        fname0 = os.path.join(work_dir, 'timeseries.h5')
        fname1 = os.path.join(work_dir, 'timeseries.h5')
        # 读取 timeseries.h5 的属性，用来判断卫星平台、是否已地理编码等信息。
        atr = readfile.read_attribute(fname0)

        phase_correction_steps = [
            'correct_LOD',
            'correct_SET',
            'correct_ionosphere',
            'correct_troposphere',
            'deramp',
            'correct_topography',
        ]

        # loop for all steps
        steps = dict()
        for sname in phase_correction_steps:
            # 默认情况下输入文件与输出文件相同，表示该步骤无需生成新文件。
            # fname0 == fname1 if no valid correction method is set.
            fname0 = fname1
            if sname == 'correct_LOD':
                # Envisat 平台才需要 LOD 改正，输出文件名增加 _LOD 后缀。
                if atr['PLATFORM'].lower().startswith('env'):
                    fname1 = f'{os.path.splitext(fname0)[0]}_LOD.h5'

            elif sname == 'correct_SET':
                method = template['mintpy.solidEarthTides']
                if method:
                    # 如果模板开启固体地球潮汐改正，输出文件名增加 _SET 后缀。
                    fname1 = f'{os.path.splitext(fname0)[0]}_SET.h5'

            elif sname == 'correct_ionosphere':
                method = template['mintpy.ionosphericDelay.method']
                if method:
                    if method == 'split_spectrum':
                        # 电离层改正输出文件名增加 _ion 后缀。
                        fname1 = f'{os.path.splitext(fname0)[0]}_ion.h5'
                    else:
                        msg = f'un-recognized ionospheric correction method: {method}'
                        raise ValueError(msg)

            elif sname == 'correct_troposphere':
                method = template['mintpy.troposphericDelay.method']
                model  = template['mintpy.troposphericDelay.weatherModel']
                if method:
                    if method == 'height_correlation':
                        # 高程相关对流层改正输出文件名增加 _tropHgt 后缀。
                        fname1 = f'{os.path.splitext(fname0)[0]}_tropHgt.h5'

                    elif method == 'gacos':
                        # GACOS 对流层改正输出文件名增加 _GACOS 后缀。
                        fname1 = f'{os.path.splitext(fname0)[0]}_GACOS.h5'

                    elif method == 'pyaps':
                        # PyAPS 会把天气模型名写入输出文件名，例如 _ERA5.h5。
                        fname1 = f'{os.path.splitext(fname0)[0]}_{model}.h5'

                    else:
                        msg = f'un-recognized tropospheric correction method: {method}'
                        raise ValueError(msg)

            elif sname == 'deramp':
                method = template['mintpy.deramp']
                if method:
                    if method in RAMP_LIST:
                        # 去相位坡度后的文件名增加 _ramp 后缀。
                        fname1 = f'{os.path.splitext(fname0)[0]}_ramp.h5'
                    else:
                        msg = f'un-recognized phase ramp type: {method}'
                        msg += f'\navailable ramp types:\n{RAMP_LIST}'
                        raise ValueError(msg)

            elif sname == 'correct_topography':
                method = template['mintpy.topographicResidual']
                if method:
                    # DEM 残差改正后的文件名增加 _demErr 后缀。
                    fname1 = f'{os.path.splitext(fname0)[0]}_demErr.h5'

            step = dict()
            # 每个步骤记录 input 和 output，后面的 run_xxx 方法会根据它们判断是否需要运行。
            step['input'] = fname0
            step['output'] = fname1
            steps[sname] = step

        # step - reference_date
        fnames = [steps[sname]['output'] for sname in phase_correction_steps]
        fnames += [steps[sname]['input'] for sname in phase_correction_steps]
        # set 去重，sorted 排序；这样 reference_date 可以处理所有可能产生过的时序文件。
        fnames = sorted(list(set(fnames)))
        step = dict()
        step['input'] = fnames
        steps['reference_date'] = step

        # velocity 和 geocode 都使用参考日期处理后的最终时序文件作为输入。
        # step - velocity / geocode
        step = dict()
        step['input'] = steps['reference_date']['input'][-1]
        steps['velocity'] = step
        steps['geocode'] = step

        # 若当前数据仍处于雷达坐标，HDF-EOS5 输出应使用 geo 目录中的地理编码结果。
        # step - hdfeos5
        if 'Y_FIRST' not in atr.keys():
            step = dict()
            fdir = os.path.dirname(steps['reference_date']['input'][-1])
            fbase = os.path.basename(steps['reference_date']['input'][-1])
            step['input'] = os.path.join(fdir, f'geo/geo_{fbase}')
        steps['hdfeos5'] = step
        return steps


    def run_local_oscillator_drift_correction(self, step_name):
        """Correct local oscillator drift (LOD).
        Automatically applied for Envisat data.
        Automatically skipped for all the other data.
        """
        # 本地振荡器漂移改正主要针对 Envisat 数据，其它卫星通常跳过。
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        # get_timeseries_filename() 返回一个字典，里面告诉当前步骤的输入文件和输出文件。
        fnames = self.get_timeseries_filename(self.template, self.workDir)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            iargs = [in_file, geom_file, '-o', out_file]
            print('\nlocal_oscilator_drift.py', ' '.join(iargs))
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                # 调用 local_oscilator_drift.py 执行 Envisat LOD 改正。
                import mintpy.cli.local_oscilator_drift
                mintpy.cli.local_oscilator_drift.main(iargs)
        else:
            # 输入和输出文件名相同，说明这个步骤不需要生成新文件。
            atr = readfile.read_attribute(in_file)
            sat = atr.get('PLATFORM', None)
            print(f'No local oscillator drift correction is needed for {sat}.')


    def run_solid_earth_tides_correction(self, step_name):
        """Correct solid Earth tides (SET)."""
        # 根据几何文件和模板设置执行固体地球潮汐改正。
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        fnames = self.get_timeseries_filename(self.template, self.workDir)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']

        if in_file != out_file:
            iargs = [in_file, '-g', geom_file, '-o', out_file, '--update']
            print('\nsolid_earth_tides.py', ' '.join(iargs))
            if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                # 调用 solid_earth_tides.py 给时序文件做固体地球潮汐改正。
                import mintpy.cli.solid_earth_tides
                mintpy.cli.solid_earth_tides.main(iargs)
        else:
            print('No solid Earth tides correction.')


    def run_ionospheric_delay_correction(self, step_name):
        """Correct ionospheric delays."""
        # 电离层改正依赖加载阶段生成的 ionosphere stack。
        iono_stack_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[3]

        fnames = self.get_timeseries_filename(self.template, self.workDir)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            method = self.template['mintpy.ionosphericDelay.method']

            # 当前工作流支持 split_spectrum 方法。
            # range split spectrum (Fattahi et al., 2017; Liang et al. 2018; 2019)
            if method == 'split_spectrum':
                print(f'ionospheric delay correction with {method} approach')
                iargs = ['-t', self.templateFile, '-f', in_file, '-o', out_file,
                         '--iono-stack-file', iono_stack_file]
                print('\niono_split_spectrum.py', ' '.join(iargs))
                # 调用 iono_split_spectrum.py 使用分裂谱方法改正电离层延迟。
                import mintpy.cli.iono_split_spectrum
                mintpy.cli.iono_split_spectrum.main(iargs)
        else:
            print('No ionospheric delay correction.')


    def run_tropospheric_delay_correction(self, step_name):
        """Correct tropospheric delays."""
        # 对流层延迟改正支持高程相关、GACOS 和 PyAPS 天气再分析三类方法。
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        # maskTempCoh.h5 是可靠像元掩膜，用来限制对流层拟合只在可信像元上进行。
        mask_file = os.path.join(self.workDir, 'maskTempCoh.h5')

        fnames = self.get_timeseries_filename(self.template, self.workDir)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            # 从模板读取对流层改正所需参数。
            poly_order  = self.template['mintpy.troposphericDelay.polyOrder']
            tropo_model = self.template['mintpy.troposphericDelay.weatherModel'].upper()
            weather_dir = self.template['mintpy.troposphericDelay.weatherDir']
            method      = self.template['mintpy.troposphericDelay.method']

            def get_dataset_size(fname):
                # 用文件属性中的 LENGTH/WIDTH 判断两个 HDF5 数据集尺寸是否一致。
                atr = readfile.read_attribute(fname)
                return (atr['LENGTH'], atr['WIDTH'])

            # Phase/Elevation Ratio (Doin et al., 2009)
            if method == 'height_correlation':
                tropo_look = self.template['mintpy.troposphericDelay.looks']
                tropo_min_cor = self.template['mintpy.troposphericDelay.minCorrelation']
                # -g 几何文件，-p 多项式阶数，-m 掩膜文件，-l looks，-t 最小相关系数阈值。
                iargs = [in_file,
                         '-g', geom_file,
                         '-p', poly_order,
                         '-m', mask_file,
                         '-o', out_file,
                         '-l', tropo_look,
                         '-t', tropo_min_cor]
                print('tropospheric delay correction with height-correlation approach')
                print('\ntropo_phase_elevation.py', ' '.join(iargs))
                if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    # 调用 tropo_phase_elevation.py 用相位-高程相关关系估计对流层延迟。
                    import mintpy.cli.tropo_phase_elevation
                    mintpy.cli.tropo_phase_elevation.main(iargs)

            # Weather re-analysis data with iterative tropospheric decomposition (GACOS)
            # Yu et al. (2018, JGR)
            elif method == 'gacos':
                # GACOS_dir 是 GACOS 数据所在目录，由模板指定。
                GACOS_dir = self.template['mintpy.troposphericDelay.gacosDir']
                iargs = ['-f', in_file, '-g', geom_file, '-o', out_file, '--dir', GACOS_dir]
                print('tropospheric delay correction with gacos approach')
                print('\ntropo_gacos.py', ' '.join(iargs))
                if ut.run_or_skip(out_file=out_file, in_file=in_file) == 'run':
                    # 调用 tropo_gacos.py 使用 GACOS 产品改正对流层延迟。
                    import mintpy.cli.tropo_gacos
                    mintpy.cli.tropo_gacos.main(iargs)

            # Weather Re-analysis Data (Jolivet et al., 2011;2014)
            elif method == 'pyaps':
                # PyAPS 需要天气模型名称、几何文件以及天气数据保存目录。
                iargs = ['-f', in_file, '--model', tropo_model, '-g', geom_file, '-w', weather_dir]
                print('Atmospheric correction using Weather Re-analysis dataset (PyAPS, Jolivet et al., 2011)')
                print('Weather Re-analysis dataset:', tropo_model)
                # 预估/缓存的天气模型延迟文件通常放在 inputs 目录下。
                tropo_file = f'./inputs/{tropo_model}.h5'

                if ut.run_or_skip(out_file=out_file, in_file=[in_file, tropo_file]) == 'run':
                    if os.path.isfile(tropo_file) and get_dataset_size(tropo_file) == get_dataset_size(in_file):
                        # 若已有尺寸匹配的对流层延迟文件，则直接与时序文件相减。
                        iargs = [in_file, tropo_file, '-o', out_file, '--force']
                        print('--------------------------------------------')
                        print(f'Use existed tropospheric delay file: {tropo_file}')
                        print('\ndiff.py', ' '.join(iargs))
                        # diff.py 这里用于从输入时序中减去已有的对流层延迟文件。
                        import mintpy.cli.diff
                        mintpy.cli.diff.main(iargs)

                    else:
                        if tropo_model in ['ERA5']:
                            # ERA5 使用 PyAPS3 接口。
                            print('\ntropo_pyaps3.py', ' '.join(iargs))
                            # tropo_pyaps3.py 会下载/读取 ERA5 并生成对流层改正后的时序。
                            import mintpy.cli.tropo_pyaps3
                            mintpy.cli.tropo_pyaps3.main(iargs)

                        elif tropo_model in ['MERRA', 'NARR']:
                            # MERRA/NARR 保留旧版 PyAPS 处理入口。
                            print('\ntropo_pyaps.py', ' '.join(iargs))
                            # legacy 表示旧版模块；这里为了兼容 MERRA/NARR 仍然调用它。
                            import mintpy.legacy.tropo_pyaps
                            mintpy.legacy.tropo_pyaps.main(iargs)

                        else:
                            raise ValueError(f'un-recognized dataset name: {tropo_model}.')

        else:
            print('No tropospheric delay correction.')


    def run_phase_deramping(self, step_name):
        """Estimate and remove phase ramp from each acquisition."""
        # 按模板指定的 ramp 类型估计并移除每期影像中的相位坡度。
        # mask_file 指定哪些像元参与坡度估计；method 指定坡度模型类型。
        mask_file = self.template['mintpy.deramp.maskFile']
        method    = self.template['mintpy.deramp']

        fnames = self.get_timeseries_filename(self.template, self.workDir)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']
        if in_file != out_file:
            print(f'Remove for each acquisition a phase ramp: {method}')
            iargs = [in_file, '-s', method, '-m', mask_file, '-o', out_file, '--update']
            print('\nremove_ramp.py', ' '.join(iargs))
            # 调用 remove_ramp.py 去除每期影像中的平面/二次等相位坡度。
            import mintpy.cli.remove_ramp
            mintpy.cli.remove_ramp.main(iargs)
        else:
            print('No phase ramp removal.')


    def run_topographic_residual_correction(self, step_name):
        """step - correct_topography
        Topographic residual (DEM error) correction (optional).
        """
        # 估计并改正 DEM 误差造成的残余地形相位。
        geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
        fnames = self.get_timeseries_filename(self.template, self.workDir)[step_name]
        in_file = fnames['input']
        out_file = fnames['output']

        if in_file != out_file:
            iargs = [in_file, '-t', self.templateFile, '-o', out_file, '--update']
            if self.template['mintpy.topographicResidual.pixelwiseGeometry']:
                # 如果模板要求逐像元几何，就额外把几何文件通过 -g 传给 dem_error.py。
                iargs += ['-g', geom_file]
            print('\ndem_error.py', ' '.join(iargs))
            # dem_error.py 估计 DEM 残差并从时序中扣除相关相位。
            import mintpy.cli.dem_error
            mintpy.cli.dem_error.main(iargs)

        else:
            print('No topographic residual correction.')


    def run_residual_phase_rms(self, step_name):
        """Noise evaluation based on the phase residual."""
        # 若残差时序文件存在，则计算 RMS 作为噪声诊断指标。
        res_file = 'timeseriesResidual.h5'
        if os.path.isfile(res_file):
            iargs = [res_file, '-t', self.templateFile]
            print('\ntimeseries_rms.py', ' '.join(iargs))
            # timeseries_rms.py 根据残差时序计算 RMS，用于评估噪声水平。
            import mintpy.cli.timeseries_rms
            mintpy.cli.timeseries_rms.main(iargs)
        else:
            print('No residual phase file found! Skip residual RMS analysis.')


    def run_reference_date(self, step_name):
        """Change reference date for all time-series files (optional)."""
        # 将所有中间/最终时序文件统一调整到模板指定的参考日期。
        if self.template['mintpy.reference.date']:
            iargs = ['-t', self.templateFile]
            in_files = self.get_timeseries_filename(self.template, self.workDir)[step_name]['input']
            for in_file in in_files:
                # reference_date.py 支持一次传入多个时序文件，这里逐个追加到参数列表。
                iargs += [in_file]
            print('\nreference_date.py', ' '.join(iargs))
            # 调用 reference_date.py 修改时序文件的参考日期。
            import mintpy.cli.reference_date
            mintpy.cli.reference_date.main(iargs)
        else:
            print('No reference date change.')


    def run_timeseries2velocity(self, step_name):
        """Estimate average velocity from displacement time-series"""
        # 从最终位移时序估计平均速度。
        ts_file = self.get_timeseries_filename(self.template, self.workDir)[step_name]['input']
        # velocity.h5 是 MintPy 默认的平均速度输出文件。
        vel_file = os.path.join(self.workDir, 'velocity.h5')

        iargs = [ts_file, '-t', self.templateFile, '-o', vel_file, '--update']
        print('\ntimeseries2velocity.py', ' '.join(iargs))
        # timeseries2velocity.py 会对位移时序拟合速度模型，生成 velocity.h5。
        import mintpy.cli.timeseries2velocity
        mintpy.cli.timeseries2velocity.main(iargs)

        # Velocity from estimated tropospheric delays
        tropo_model = self.template['mintpy.troposphericDelay.weatherModel'].upper()
        tropo_file = os.path.join(self.workDir, f'inputs/{tropo_model}.h5')
        if os.path.isfile(tropo_file):
            # 如果存在估计的对流层延迟文件，也计算其速度趋势，便于评估大气改正影响。
            suffix = os.path.splitext(os.path.basename(tropo_file))[0]
            tropo_vel_file = f'{os.path.splitext(vel_file)[0]}{suffix}.h5'
            tropo_vel_file = os.path.join(self.workDir, tropo_vel_file)

            iargs = [tropo_file, '-t', self.templateFile, '-o', tropo_vel_file, '--update']
            # add reference info for a meaningful velocity to assess the impact of tropo delay on velocity
            # 为了让对流层延迟速度和主速度文件可比较，需要使用相同参考日期和参考像元。
            atr = readfile.read_attribute(vel_file)
            iargs += ['--ref-date', atr['REF_DATE'], '--ref-yx', atr['REF_Y'], atr['REF_X']]
            print('\ntimeseries2velocity.py', ' '.join(iargs))
            mintpy.cli.timeseries2velocity.main(iargs)


    def run_geocode(self, step_name):
        """geocode data files in radar coordinates into ./geo folder."""
        # 将雷达坐标产品地理编码到 geo 目录；已是地理坐标的数据会跳过。
        if self.template['mintpy.geocode']:
            ts_file = self.get_timeseries_filename(self.template, self.workDir)[step_name]['input']
            atr = readfile.read_attribute(ts_file)
            # Y_FIRST 是 MintPy 常用的地理坐标属性；没有它通常表示数据还在雷达坐标。
            if 'Y_FIRST' not in atr.keys():
                # 1. geocode
                out_dir = os.path.join(self.workDir, 'geo')
                os.makedirs(out_dir, exist_ok=True)

                # 对关键结果文件统一地理编码，保持输出产品坐标一致。
                geom_file, lookup_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1:3]
                in_files = [geom_file, 'temporalCoherence.h5', 'avgSpatialCoh.h5', ts_file, 'velocity.h5']
                # geocode.py 用 lookup_file 把雷达坐标文件转换到地理坐标，并输出到 geo 目录。
                iargs = in_files + ['-l', lookup_file, '-t', self.templateFile, '--outdir', out_dir, '--update']
                print('\ngeocode.py', ' '.join(iargs))
                import mintpy.cli.geocode
                mintpy.cli.geocode.main(iargs)

                # 2. generate reliable pixel mask in geo coordinate
                geom_file = os.path.join(out_dir, f'geo_{os.path.basename(geom_file)}')
                tcoh_file = os.path.join(out_dir, 'geo_temporalCoherence.h5')
                mask_file = os.path.join(out_dir, 'geo_maskTempCoh.h5')
                tcoh_min = self.template['mintpy.networkInversion.minTempCoh']

                iargs = [tcoh_file, '-m', tcoh_min, '-o', mask_file]
                # exclude pixels in shadow if shadowMask dataset is available
                if (self.template['mintpy.networkInversion.shadowMask'] is True
                        and 'shadowMask' in readfile.get_dataset_list(geom_file)):
                    # 在地理坐标下重新生成可靠像元掩膜，同时排除阴影区域。
                    iargs += ['--base', geom_file, '--base-dataset', 'shadowMask', '--base-value', '1']
                print('\ngenerate_mask.py', ' '.join(iargs))

                if ut.run_or_skip(out_file=mask_file, in_file=tcoh_file) == 'run':
                    # 调用 generate_mask.py 生成 geo_maskTempCoh.h5。
                    import mintpy.cli.generate_mask
                    mintpy.cli.generate_mask.main(iargs)

            else:
                print('dataset is geocoded, skip geocoding and continue.')
        else:
            print('geocoding is OFF')


    def run_save2google_earth(self, step_name):
        """Save velocity file in geo coordinates into Google Earth raster image."""
        # 按模板开关将速度结果保存为 Google Earth 可读的 KMZ 文件。
        if self.template['mintpy.save.kmz'] is True:
            print('creating Google Earth KMZ file for geocoded velocity file: ...')
            # input
            vel_file = os.path.join(self.workDir, 'velocity.h5')
            atr = readfile.read_attribute(vel_file)
            if 'Y_FIRST' not in atr.keys():
                # 雷达坐标速度文件需要改用地理编码后的版本。
                vel_file = os.path.join(self.workDir, 'geo/geo_velocity.h5')

            # output
            # os.path.splitext(vel_file)[0] 去掉 .h5 后缀，再加 .kmz 得到输出文件名。
            kmz_file = f'{os.path.splitext(vel_file)[0]}.kmz'
            iargs = [vel_file, '-o', kmz_file]
            print('\nsave_kmz.py', ' '.join(iargs))

            # update mode
            fbase = os.path.basename(kmz_file)
            # KMZ 文件可能已经在当前目录、geo 目录或 pic 目录中，这里都检查一遍。
            kmz_files = [i for i in [fbase, f'./geo/{fbase}', f'./pic/{fbase}']
                         if os.path.isfile(i)]
            kmz_file = kmz_files[0] if len(kmz_files) > 0 else None

            if ut.run_or_skip(out_file=kmz_file, in_file=vel_file, readable=False) == 'run':
                # save_kmz.py 把地理坐标速度文件转换为 Google Earth 可打开的 KMZ 栅格图。
                import mintpy.cli.save_kmz
                mintpy.cli.save_kmz.main(iargs)

        else:
            print('save velocity to Google Earth format is OFF.')


    def run_save2hdfeos5(self, step_name):
        """Save displacement time-series and its aux data in geo coordinate into HDF-EOS5 format"""
        # 按模板开关将时序、相干性、掩膜和几何信息打包为 HDF-EOS5 产品。
        if self.template['mintpy.save.hdfEos5'] is True:
            # input
            ts_file = self.get_timeseries_filename(self.template, self.workDir)[step_name]['input']
            # Add attributes from custom template to timeseries file
            if self.customTemplate is not None:
                # 输出前补充用户模板中的属性，方便产品元数据完整记录。
                ut.add_attribute(ts_file, self.customTemplate)

            tcoh_file = os.path.join(self.workDir, 'temporalCoherence.h5')
            scoh_file = os.path.join(self.workDir, 'avgSpatialCoh.h5')
            mask_file = os.path.join(self.workDir, 'maskTempCoh.h5')
            geom_file = ut.check_loaded_dataset(self.workDir, print_msg=False)[1]
            if 'geo' in ts_file:
                # 如果输入时序已经位于 geo 目录，同步使用 geo 目录中的辅助文件。
                tcoh_file = os.path.join(self.workDir, 'geo/geo_temporalCoherence.h5')
                scoh_file = os.path.join(self.workDir, 'geo/geo_avgSpatialCoh.h5')
                mask_file = os.path.join(self.workDir, 'geo/geo_maskTempCoh.h5')
                geom_file = os.path.join(self.workDir, f'geo/geo_{os.path.basename(geom_file)}')

            # cmd
            print('--------------------------------------------')
            # save_hdfeos5.py 需要同时知道时序、时间相干性、空间相干性、掩膜、几何和模板文件。
            iargs = [ts_file,
                     '--tc', tcoh_file,
                     '--asc', scoh_file,
                     '-m', mask_file,
                     '-g', geom_file,
                     '-t', self.templateFile]
            print('\nsave_hdfeos5.py', ' '.join(iargs))

            # output (check existing file)
            atr = readfile.read_attribute(ts_file)
            # sensor.get_unavco_mission_name() 根据文件属性推断 UNAVCO/HDF-EOS5 产品所需的卫星任务名。
            SAT = sensor.get_unavco_mission_name(atr)
            # ut.get_file_list('SAT_*.he5') 查找是否已经有对应卫星名前缀的 HDF-EOS5 输出文件。
            hdfeos5_files = ut.get_file_list(f'{SAT}_*.he5')
            hdfeos5_file = hdfeos5_files[0] if len(hdfeos5_files) > 0 else None

            if ut.run_or_skip(out_file=hdfeos5_file,
                              in_file=[ts_file, tcoh_file,
                                       scoh_file, mask_file,
                                       geom_file]) == 'run':
                # save_hdfeos5.py 生成标准 HDF-EOS5 格式产品。
                import mintpy.cli.save_hdfeos5
                mintpy.cli.save_hdfeos5.main(iargs)

        else:
            print('save time-series to HDF-EOS5 format is OFF.')


    def run(self, steps):
        """run the chosen steps."""
        # 按解析出的步骤列表顺序执行，每个步骤映射到对应的成员函数。
        for sname in steps:
            print(f'\n\n******************** step - {sname} ********************')

            # 这一长串 if/elif 是“调度表”：根据步骤名 sname，调用对应的类方法。
            # 例如 sname == 'load_data' 时，就调用 self.run_load_data(sname)。
            if sname == 'load_data':
                self.run_load_data(sname)

            elif sname == 'modify_network':
                self.run_network_modification(sname)

            elif sname == 'reference_point':
                self.run_reference_point(sname)

            elif sname == 'quick_overview':
                self.run_quick_overview(sname)

            elif sname == 'correct_unwrap_error':
                self.run_unwrap_error_correction(sname)

            elif sname == 'invert_network':
                self.run_network_inversion(sname)

            elif sname == 'correct_LOD':
                self.run_local_oscillator_drift_correction(sname)

            elif sname == 'correct_SET':
                self.run_solid_earth_tides_correction(sname)

            elif sname == 'correct_ionosphere':
                self.run_ionospheric_delay_correction(sname)

            elif sname == 'correct_troposphere':
                self.run_tropospheric_delay_correction(sname)

            elif sname == 'deramp':
                self.run_phase_deramping(sname)

            elif sname == 'correct_topography':
                self.run_topographic_residual_correction(sname)

            elif sname == 'residual_RMS':
                self.run_residual_phase_rms(sname)

            elif sname == 'reference_date':
                self.run_reference_date(sname)

            elif sname == 'velocity':
                self.run_timeseries2velocity(sname)

            elif sname == 'geocode':
                self.run_geocode(sname)

            elif sname == 'google_earth':
                self.run_save2google_earth(sname)

            elif sname == 'hdfeos5':
                self.run_save2hdfeos5(sname)


    def plot_result(self, print_aux=True):
        """Plot data files and save to figures in pic folder"""

        # 将主要中间结果和最终结果批量绘图，并统一移动到 pic 目录。
        print('\n******************** plot & save to pic ********************')

        # 从模板读取绘图和计算资源相关参数。
        tropo_model = self.template['mintpy.troposphericDelay.weatherModel'].upper()
        max_plot_memory = abs(float(self.template['mintpy.plot.maxMemory']))
        max_memory = abs(float(self.template['mintpy.compute.maxMemory']))
        fig_dpi = int(self.template['mintpy.plot.dpi'])

        # 再次检查核心输入文件位置，绘图时会用到它们。
        stack_file, geom_file, lookup_file, ion_file = ut.check_loaded_dataset(
            self.workDir,
            print_msg=False)[:4]
        mask_file = os.path.join(self.workDir, 'maskTempCoh.h5')
        geo_dir = os.path.join(self.workDir, 'geo')
        pic_dir = os.path.join(self.workDir, 'pic')

        # 打印 view 命令时使用相对路径，输出更短也更容易阅读。
        # use relative path for shorter and cleaner printout view command
        stack_file  = os.path.relpath(stack_file)  if stack_file  else stack_file
        ion_file    = os.path.relpath(ion_file)    if ion_file    else ion_file
        geom_file   = os.path.relpath(geom_file)   if geom_file   else geom_file
        lookup_file = os.path.relpath(lookup_file) if lookup_file else lookup_file
        mask_file   = os.path.relpath(mask_file)   if mask_file   else mask_file
        geo_dir     = os.path.relpath(geo_dir)     if geo_dir     else geo_dir

        # view options
        # for each element list:
        # the 1st item is the data file
        # the 2nd item is the dataset if applicable
        # opt4ts 是绘制时序文件时常用的一组 view.py 参数。
        opt4ts = ['--noaxis', '-u', 'cm', '--wrap', '--wrap-range', '-5', '5']
        # iargs_list0 中每个小列表都代表一次 view.py 调用的参数。
        # 例如 ['velocity.h5', '--dem', geom_file, '--mask', mask_file] 表示绘制 velocity.h5。
        iargs_list0 = [
            # key files
            ['velocity.h5',          '--dem', geom_file, '--mask', mask_file],
            ['temporalCoherence.h5', '-c', 'gray', '-v', '0', '1'],
            ['maskTempCoh.h5',       '-c', 'gray', '-v', '0', '1'],

            # geometry
            [geom_file],
            [lookup_file],

            # ifgramStack
            [stack_file, 'unwrapPhase-',      '--zero-mask', '--wrap', '-c', 'cmy'],
            [stack_file, 'unwrapPhase-',      '--zero-mask'],
            [stack_file, 'coherence-',        '--mask', 'no', '-v', '0', '1'],
            [stack_file, 'connectComponent-', '--mask', 'no'],

            # ifgramStack - unwrapping error correction
            [stack_file, 'unwrapPhase_bridging-',              '--zero-mask'],
            [stack_file, 'unwrapPhase_phaseClosure-',          '--zero-mask'],
            [stack_file, 'unwrapPhase_bridging_phaseClosure-', '--zero-mask'],

            # ifgramStack - auxliary files
            ['avgPhaseVelocity.h5'] ,
            ['avgSpatialCoh.h5', '-c', 'gray', '-v', '0', '1'],
            ['maskConnComp.h5',  '-c', 'gray', '-v', '0', '1'],

            # time-series
            ['timeseries.h5'] + opt4ts,
            ['timeseries_*.h5'] + opt4ts,

            # files from geocoding
            [os.path.join(geo_dir, 'geo_maskTempCoh.h5'),       '-c', 'gray'],
            [os.path.join(geo_dir, 'geo_temporalCoherence.h5'), '-c', 'gray'],
            [os.path.join(geo_dir, 'geo_avgSpatialCoh.h5'),     '-c', 'gray'],
            [os.path.join(geo_dir, 'geo_velocity.h5'),          'velocity'],
            [os.path.join(geo_dir, 'geo_timeseries*.h5')] + opt4ts,

            # all the other files
            [f'velocity{tropo_model}.h5', '--mask', 'no'],
            ['numInvIfgram.h5',           '--mask', 'no'],
        ]

        if ion_file:
            # 如果存在电离层栈，也把相关数据集加入绘图列表。
            iargs_list0 += [
                [ion_file, 'unwrapPhase-', '--zero-mask', '--wrap', '-c', 'cmy'],
                [ion_file, 'unwrapPhase-', '--zero-mask'],
                [ion_file, 'coherence-',   '--mask', 'no', '-v', '0', '1'],
            ]

        # translate element list whose file path has *
        iargs_list = []
        for iargs in iargs_list0:
            fname, args = iargs[0], iargs[1:]
            if not fname:
                # 如果文件名是 None 或空字符串，说明这个输入不存在，直接跳过。
                continue

            if '*' in fname:
                # 展开通配符路径，逐个生成 view.py 参数列表。
                # glob.glob() 会把 'timeseries_*.h5' 这种通配符展开成真实文件名列表。
                fnames = sorted(glob.glob(fname))
                if len(fnames) > 0:
                    for fname in fnames:
                        iargs_list.append([fname] + args)

            elif iargs not in iargs_list:
                iargs_list.append(iargs)

        # remove element list - file does not exists
        # 只保留真实存在的文件，避免 view.py 因找不到文件而报错。
        iargs_list = [iargs for iargs in iargs_list if os.path.isfile(iargs[0])]

        # remote element list - file is ifgramStack and the dataset of interest does not exist
        # 对 ifgramStack.h5 这种一个文件包含多个数据集的情况，检查要绘制的数据集是否真的存在。
        stack_dset_list = readfile.get_dataset_list(stack_file)
        stack_dset_list += [f'{x}-' for x in stack_dset_list]
        iargs_list = [iargs for iargs in iargs_list
                      if (iargs[0] != stack_file
                          or (iargs[0] == stack_file
                              and iargs[1] in stack_dset_list))]

        # add the following common options to all element lists
        # opt_common 是所有 view.py 调用共用的绘图选项，例如 dpi、无界面显示、内存限制等。
        opt_common = ['--dpi', str(fig_dpi), '--noverbose', '--nodisplay',
                      '--update', '--memory', str(max_plot_memory)]
        # 把通用参数加到每一次 view.py 调用前面。
        iargs_list = [opt_common + iargs for iargs in iargs_list]

        # run view
        start_time = time.time()
        run_parallel = False
        cluster_type = self.template['mintpy.compute.cluster']
        if cluster_type:
            # 若模板指定并行计算，则根据 CPU 数和内存限制确定实际并行度。
            num_workers = self.template['mintpy.compute.numWorker']
            # DaskCluster.format_num_worker() 会把模板里的 worker 数量配置转换成可用的数字。
            num_workers = cluster.DaskCluster.format_num_worker(cluster_type, num_workers)

            # limit number of parallel processes based on available CPU
            num_cores, run_parallel, Parallel, delayed = ut.check_parallel(
                len(iargs_list),
                print_msg=False,
                maxParallelNum=num_workers)

            # limit number of parallel processes based on memory limit
            # default view.py call could use up to 1.5 GB reserved memory (~700 MB actual memory)
            # scale based on utils.plot.auto_multilook_num()
            plot_memory = 1.5 if 2.0 < max_plot_memory <= 4.0 else 1.5 * (max_plot_memory / 4.0)
            num_cores = min(num_cores, max(int(max_memory / plot_memory), 1))

        import mintpy.cli.view
        if run_parallel and num_cores > 1:
            print(f"parallel processing using {num_cores} cores ...")
            # Parallel/delayed 来自 joblib；这里把多次 view.py 绘图任务并行运行。
            Parallel(n_jobs=num_cores)(delayed(mintpy.cli.view.main)(iargs) for iargs in iargs_list)
        else:
            for iargs in iargs_list:
                # 串行模式：一个一个调用 view.py 绘图。
                mintpy.cli.view.main(iargs)

        # copy text files to pic
        print('copy *.txt files into ./pic directory.')
        # glob.glob('*.txt') 查找当前目录下所有 txt 文件。
        tfiles = glob.glob('*.txt')
        for tfile in tfiles:
            # 把文本诊断文件复制到 pic 目录，方便集中查看。
            shutil.copy2(tfile, pic_dir)

        # move picture files to pic
        print('move *.png/pdf/kmz files to ./pic directory.')
        pfiles  = glob.glob('*.png')
        pfiles += glob.glob('*.pdf')
        pfiles += glob.glob('*.kmz')
        pfiles += glob.glob(os.path.join(geo_dir, '*.kmz'))
        for pfile in pfiles:
            # shutil.move() 移动文件到 pic 目录；os.path.basename() 保留原文件名。
            shutil.move(pfile, os.path.join(pic_dir, os.path.basename(pfile)))

        # time info
        # time.time() 返回当前时间戳；与 start_time 相减得到绘图耗时秒数。
        # divmod(seconds, 60) 把秒数拆成分钟和剩余秒数。
        m, s = divmod(time.time()-start_time, 60)
        print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

        # message for more visualization scripts
        msg = """Explore more info & visualization options with the following scripts:
        info.py                    # check HDF5 file structure and metadata
        view.py                    # 2D map view
        tsview.py                  # 1D point time-series (interactive)
        transect.py                # 1D profile (interactive)
        plot_coherence_matrix.py   # plot coherence matrix for one pixel (interactive)
        plot_network.py            # plot network configuration of the dataset
        plot_transection.py        # plot 1D profile along a line of a 2D matrix (interactive)
        save_kmz.py                # generate Google Earth KMZ file in raster image
        save_kmz_timeseries.py     # generate Google Earth KMZ file in points for time-series (interactive)
        """
        if print_aux:
            print(msg)


    def close(self, normal_end=True):
        # 流程结束后回到启动程序时所在的目录。
        # go back to original directory
        print('Go back to directory:', self.cwd)
        # os.chdir(self.cwd) 切回最开始的目录，避免程序结束后还停留在工作目录中。
        os.chdir(self.cwd)
        # message
        if normal_end:
            msg  = '\n################################################'
            msg += '\n   Normal end of smallbaselineApp processing!'
            msg += '\n################################################'
            print(msg)


def run_smallbaselineApp(inps):
    """Run the small baseline time series analsysis workflow."""
    # smallbaselineApp 的程序化入口：打开工作流、执行步骤、按需绘图并收尾。
    start_time = time.time()

    # open and run
    # 创建 TimeSeriesAnalysis 对象；这会自动调用 __init__ 保存模板路径和工作目录。
    app = TimeSeriesAnalysis(inps.customTemplateFile, inps.workDir)
    # app.open() 负责进入工作目录、准备模板、读取配置。
    app.open()
    # app.run(...) 根据 CLI 文件解析出来的 runSteps 逐步执行处理流程。
    app.run(steps=inps.runSteps)

    # plot if:
    # a) --plot in command line OR
    # b) template['mintpy.plot'] = yes AND runSteps > 1
    # 如果用户明确要求绘图，或模板开启绘图且本次运行了多个步骤，就调用 plot_result()。
    if inps.plot or (app.template['mintpy.plot'] and len(inps.runSteps) > 1):
        app.plot_result()

    # close
    # close() 会切回原目录，并打印正常结束信息。
    app.close()

    # used time
    # 统计整个 smallbaselineApp 工作流从开始到结束的总耗时。
    m, s = divmod(time.time()-start_time, 60)
    print(f'Time used: {m:02.0f} mins {s:02.1f} secs\n')

    return
