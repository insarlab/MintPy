#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################

# os 用来处理文件路径，例如取当前脚本文件名、判断模板文件名。
import os
# sys 用来读取命令行参数和退出程序；sys.exit(0) 表示正常退出，sys.exit(1) 表示出错退出。
import sys

# auto_path 保存不同处理软件（ISCE/GAMMA 等）的默认输入路径规则。
from mintpy.defaults import auto_path
# get_template_content() 会读取 MintPy 内置的模板说明文本，显示在命令行帮助中。
from mintpy.defaults.template import get_template_content
# create_argument_parser() 是 MintPy 封装的命令行解析器创建函数。
from mintpy.utils.arg_utils import create_argument_parser

#################################################################
# DEFAULT_TEMPLATE 是 `load_data.py -H` 时显示的示例模板。
# """...""" 是多行字符串；format(...) 把不同处理器的自动路径说明拼进去。
DEFAULT_TEMPLATE = """template:
########## 1. Load Data (--load to exit after this step)
{}\n
{}\n
{}""".format(
    auto_path.AUTO_PATH_GAMMA,
    auto_path.AUTO_PATH_ISCE_STRIPMAP,
    auto_path.AUTO_PATH_ISCE_TOPS,
)

# TEMPLATE 是 load_data 对应的完整模板配置说明。
TEMPLATE = get_template_content('load_data')

# NOTE 是补充说明，提醒用户哪些数据集是必需的、文件名里需要包含什么日期信息。
NOTE = """NOTE:
  For interferogram, unwrapPhase is required, the other dataset are optional, including coherence, connectComponent, wrapPhase, etc.
  The unwrapPhase metadata file requires DATE12 attribute in YYMMDD-YYMMDD format.
  All path of data file must contain the reference and secondary date, either in file name or folder name.
"""

# EXAMPLE 是命令行示例，会显示在 `load_data.py -h` 的帮助信息里。
EXAMPLE = """example:
  # MUST run in the mintpy working directory!

  # show example template file for ISCE/ROI_PAC/GAMMA products
  load_data.py -H

  # load & write the following HDF5 files:
  # ./inputs/ifgramStack.h5   for interferogram        stack
  # ./inputs/ionStack.h5      for ionosphere           stack
  # ./inputs/offsetStack.h5   for range/azimuth offset stack
  # ./inputs/geometryRadar.h5 for geometry in radar coordinates
  # ./inputs/geometryGeo.h5   for geometry in geo   coordinates
  load_data.py -t smallbaselineApp.cfg
  load_data.py -t smallbaselineApp.cfg GalapagosSenDT128.txt --project GalapagosSenDT128

  # load geometry ONLY
  smallbaselineApp.py SaltonSeaSenDT173.txt -g
  load_data.py -t smallbaselineApp.cfg --geom
"""


def create_parser(subparsers=None):
    """Create command line parser."""
    # 创建命令行解析器：把终端输入的 `load_data.py -t xxx.cfg` 转换成 Python 变量。
    synopsis = 'Load stacks of interferograms to HDF5 files'
    epilog = TEMPLATE + '\n' + NOTE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # extra help
    # action='store_true' 表示用户输入 -H 时，inps.print_example_template 变成 True。
    parser.add_argument('-H', dest='print_example_template', action='store_true',
                        help='Print/Show the example template file for loading.')

    # input files
    # nargs='+' 表示 -t 后面可以跟一个或多个模板文件。
    parser.add_argument('-t', '--template', dest='template_file', type=str, nargs='+',
                        help='template file(s) with path info.')
    # --geom 只加载几何文件，不加载干涉图栈；常用于先检查几何输入是否正确。
    parser.add_argument('--geom','--geometry', dest='only_load_geometry', action='store_true',
                        help='Load the geometry file(s) ONLY.')

    # options from template file name & content
    # PROJECT_NAME 会写入输出 HDF5 元数据，便于后续产品识别。
    parser.add_argument('--project', type=str, dest='PROJECT_NAME',
                        help='project name of dataset for INSARMAPS Web Viewer')
    # --enforce 会把 updateMode 设为 False，意思是强制重新加载，不用已有文件跳过。
    parser.add_argument('--enforce', '-f', dest='updateMode', action='store_false',
                        help='Disable the update mode, or skip checking dataset already loaded.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    # parse_args() 解析参数；iargs 为 None 时读取真实命令行，传列表时方便其它代码或测试调用。
    inps = parser.parse_args(args=iargs)

    # check: -H option
    if inps.print_example_template:
        # 只打印示例模板并退出，不执行真正的数据加载。
        print(DEFAULT_TEMPLATE)
        sys.exit(0)

    # check: -t/--template option
    # -t option is required AND
    # smallbaselineApp.cfg file is required
    if not inps.template_file:
        # 如果没有 -t/--template，就打印用法并退出。
        parser.print_usage()
        script_name = os.path.basename(__file__)
        print(f'{script_name}: error: -t/--template option is required.')
        print(f'run {script_name} -H to show the example template file.')
        sys.exit(1)

    elif all(not x.endswith('smallbaselineApp.cfg') for x in inps.template_file):
        # load_data 至少需要默认模板 smallbaselineApp.cfg，因为很多默认路径和选项都来自它。
        script_name = os.path.basename(__file__)
        print(f'{script_name}: error: at least smallbaselineApp.cfg file is required for -t/--template option.')
        sys.exit(1)

    # 返回解析后的参数对象，后面的 main() 会把它交给核心 load_data() 函数。
    return inps



#################################################################
def main(iargs=None):
    # parse
    # main() 是 CLI 入口函数：先解析命令行参数。
    inps = cmd_line_parse(iargs)

    # import
    # 延迟导入核心函数，只有真正要运行时才加载 src/mintpy/load_data.py。
    from mintpy.load_data import load_data

    # run
    # 调用核心模块的 load_data(inps)，把外部数据读入 MintPy 标准 HDF5 文件。
    load_data(inps)


#################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
