#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


# datetime 是 Python 标准库，专门用来处理日期和时间；这里用它打印程序开始运行的时间。
import datetime
# os 是 Python 标准库，专门用来和操作系统交互，比如判断文件是否存在、拼接路径、获取当前目录。
import os
# sys 是 Python 标准库，专门用来访问 Python 解释器相关信息；这里用 sys.argv 读取命令行参数，用 sys.exit 退出程序。
import sys

# mintpy 是本项目的主包；这里用它读取软件版本、默认配置文件所在位置等信息。
import mintpy
# STEP_LIST 是 MintPy 默认模板中定义好的处理步骤顺序，例如 load_data、modify_network、velocity 等。
from mintpy.defaults.template import STEP_LIST
# create_argument_parser 是 MintPy 自己封装的 argparse 工具，用来创建风格统一的命令行参数解析器。
from mintpy.utils.arg_utils import create_argument_parser

##########################################################################
# 命令行帮助文本：说明可通过 --start/--end/--dostep 选择运行的处理步骤。
# """...""" 是 Python 的多行字符串，适合保存较长的帮助说明。
# .format(...) 会把 STEP_LIST 分成几段填入上面的 {}，这样 help 页面里能显示所有步骤名称。
STEP_HELP = """Command line options for steps processing with names chosen from the following list:

{}
{}
{}
{}

In order to use either --start or --dostep, it is necessary that a
previous run was done using one of the steps options to process at least
through the step immediately preceding the starting step of the current run.
""".format(STEP_LIST[0:5], STEP_LIST[5:10], STEP_LIST[10:16], STEP_LIST[16:])

REFERENCE = """reference:
  Yunjun, Z., Fattahi, H., and Amelung, F. (2019), Small baseline InSAR time series analysis:
    Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
    doi:10.1016/j.cageo.2019.104331.
"""

# EXAMPLE 保存命令行使用示例，会显示在 -h/--help 帮助信息的末尾。
EXAMPLE = """example:
  smallbaselineApp.py                         # run with default template 'smallbaselineApp.cfg'
  smallbaselineApp.py <custom_template>       # run with default and custom templates
  smallbaselineApp.py -h / --help             # help
  smallbaselineApp.py -H                      # print    default template options
  smallbaselineApp.py -g                      # generate default template if it does not exist
  smallbaselineApp.py -g <custom_template>    # generate/update default template based on custom template
  smallbaselineApp.py --plot                  # plot results w/o run [to populate the 'pic' folder after failed runs]

  # step processing with --start/stop/dostep options
  smallbaselineApp.py GalapagosSenDT128.txt --dostep velocity  # run step 'velocity' only
  smallbaselineApp.py GalapagosSenDT128.txt --end load_data    # end run after step 'load_data'
"""


def create_parser(subparsers=None):
    # 创建 smallbaselineApp.py 的命令行参数解析器。
    # “解析器”可以理解为：把用户在终端输入的文字参数，转换成 Python 里可访问的变量。
    synopsis = 'Routine Time Series Analysis for Small Baseline InSAR Stack'
    # epilog 是帮助信息最后显示的内容，这里把参考文献和示例命令拼在一起。
    epilog = REFERENCE + '\n' + EXAMPLE
    # __name__ 是 Python 自动提供的模块名；split('.')[-1] 取最后一段，作为命令名称。
    name = __name__.split('.')[-1]
    # 调用 MintPy 的 create_argument_parser() 方法，得到一个 parser 对象。
    # 后面多次调用 parser.add_argument()，就是不断告诉它“这个命令支持哪些参数”。
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # 位置参数：用户自定义模板文件；可选参数：工作目录。
    # nargs='?' 表示这个位置参数可以不填；如果用户不填，就使用当前目录里的 smallbaselineApp.cfg。
    parser.add_argument('customTemplateFile', nargs='?',
                        help='custom template with option settings.\n' +
                             "ignored if the default smallbaselineApp.cfg is input.")
    # --dir/--work-dir 是同一个参数的两个名字；dest='workDir' 表示解析后保存到 inps.workDir。
    parser.add_argument('--dir', '--work-dir', dest='workDir', default='./',
                        help='work directory, (default: %(default)s).')

    # 辅助选项：生成模板、打印模板、查看版本或仅绘图。
    # action='store_true' 表示用户写了这个选项时，对应变量就是 True；没写时就是 False。
    parser.add_argument('-g', dest='generate_template', action='store_true',
                        help='generate default template (if it does not exist) and exit.')
    parser.add_argument('-H', dest='print_template', action='store_true',
                        help='print the default template file and exit.')
    parser.add_argument('-v','--version', action='store_true', help='print software version and exit')

    parser.add_argument('--plot', dest='plot', action='store_true',
                        help='plot results [only] without running smallbaselineApp.')

    # 步骤控制选项：指定起止步骤，或只运行某一个步骤。
    # add_argument_group() 会在 help 信息中单独分组显示这几个步骤相关参数。
    step = parser.add_argument_group('steps processing (start/end/dostep)', STEP_HELP)
    # default=STEP_LIST[0] 表示默认从第一个步骤开始运行。
    step.add_argument('--start', dest='startStep', metavar='STEP', default=STEP_LIST[0],
                      help='start processing at the named step (default: %(default)s).')
    # default=STEP_LIST[-1] 表示默认运行到最后一个步骤；[-1] 是 Python 中“最后一个元素”的写法。
    step.add_argument('--end','--stop', dest='endStep', metavar='STEP',  default=STEP_LIST[-1],
                      help='end processing at the named step (default: %(default)s)')
    # --dostep 用来只运行一个步骤，后面会把 startStep 和 endStep 都改成这个步骤。
    step.add_argument('--dostep', dest='doStep', metavar='STEP',
                      help='run processing at the named step only')

    # 返回 parser，供 cmd_line_parse() 使用。
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # 解析命令行输入，返回后续工作流需要的参数对象。
    # parse
    parser = create_parser()
    # parser.parse_args() 会真正解析命令行参数。
    # 如果 iargs 是 None，表示从真实终端命令读取；如果传入列表，常用于测试或其它 Python 代码调用。
    inps = parser.parse_args(args=iargs)

    # 保存原始参数列表，用于判断用户是否只指定了特殊选项。
    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    # sys.argv 是完整命令行参数列表，第 0 个通常是脚本名，所以用 [1:] 去掉脚本名。
    inps.argv = iargs or sys.argv[1:]

    # check
    # os.path.dirname(mintpy.__file__) 得到 mintpy 包所在目录。
    # os.path.join(...) 用正确的路径分隔符拼接路径，避免手写 '/' 在不同系统上出错。
    template_file = os.path.join(os.path.dirname(mintpy.__file__), 'defaults/smallbaselineApp.cfg')

    # -H 只打印默认模板内容并退出，不进入处理流程。
    # check: -H option (print default template)
    if inps.print_template:
        # with open(...) as f 是安全打开文件的写法；代码块结束后文件会自动关闭。
        with open(template_file) as f:
            lines = f.read()
        try:
            # 如果安装了 rich，则用语法高亮显示 cfg 文件。
            # syntax highlight via rich
            from rich.console import Console
            from rich.syntax import Syntax
            # Console 是 rich 的输出对象；Syntax 负责按 cfg 语法给文本上色。
            console = Console()
            console.print(Syntax(lines, "cfg", background_color='default'))
        except ImportError:
            # 如果没有安装 rich，就退回普通 print，保证功能不依赖额外包。
            print(lines)
        # sys.exit(0) 表示正常退出程序，不再继续执行后面的处理流程。
        sys.exit(0)

    # -v 只打印 MintPy 版本信息并退出。
    # check: -v option (print software version)
    if inps.version:
        # mintpy.version.version_description 是 MintPy 的版本说明字符串。
        print(mintpy.version.version_description)
        sys.exit(0)

    # check: existence of input template files
    if (not inps.customTemplateFile
            and not os.path.isfile(os.path.basename(template_file))
            and not inps.generate_template):
        # 没有用户模板、当前目录也没有默认模板，并且不是生成模板模式时，无法继续运行。
        # os.path.basename(template_file) 只取文件名 smallbaselineApp.cfg，不带前面的目录。
        # os.path.isfile(...) 用来判断这个文件在当前目录下是否真实存在。
        parser.print_usage()
        print(EXAMPLE)

        msg = "no template file found! It requires:"
        msg += "\n  1) input a custom template file, OR"
        msg += "\n  2) there is a default template 'smallbaselineApp.cfg' in current directory."
        # raise SystemExit 会终止程序，并把错误信息显示给用户。
        raise SystemExit(f'ERROR: {msg}')

    # check: custom input template file
    if inps.customTemplateFile:
        # 自定义模板必须存在，并转换为绝对路径，避免后续切换工作目录后路径失效。
        # check the existence
        if not os.path.isfile(inps.customTemplateFile):
            # FileNotFoundError 是 Python 内置异常，用来表示文件不存在。
            raise FileNotFoundError(inps.customTemplateFile)

        # os.path.abspath(...) 把相对路径转换成绝对路径，例如 a.txt 变成 /home/.../a.txt。
        inps.customTemplateFile = os.path.abspath(inps.customTemplateFile)

        # 如果用户传入的其实是默认模板，则按无自定义模板处理。
        # ignore if it is the default template file (smallbaselineApp.cfg)
        # os.path.basename(...) 只比较文件名，不比较目录。
        if os.path.basename(inps.customTemplateFile) == os.path.basename(template_file):
            inps.customTemplateFile = None

    # 单独使用 --plot 时只绘图，不运行任何处理步骤。
    # check: --plot option
    if inps.argv == ['--plot']:
        plot_only = True
        print('plot smallbaselineApp results without run.')
    else:
        plot_only = False

    # 根据 start/end/dostep 选项计算本次需要执行的步骤列表。
    # check: --start/end/dostep options
    # read_inps2run_steps() 是本文件下面定义的函数，负责把用户输入转换成真正要运行的步骤列表。
    inps.runSteps = read_inps2run_steps(inps, step_list=STEP_LIST, plot_only=plot_only)

    # 返回 inps；后面的 main() 会把它交给真正的处理工作流。
    return inps


def read_inps2run_steps(inps, step_list, plot_only=False):
    """read/get run_steps from input arguments."""
    # 将命令行指定的步骤范围转换成有序步骤列表。
    # check: if start/end/do step input is valid
    for key in ['startStep', 'endStep', 'doStep']:
        # vars(inps) 会把 argparse 解析出来的对象转换成字典，便于用字符串 key 取值。
        value = vars(inps)[key]
        if value and value not in step_list:
            # 用户输入的步骤名必须来自默认模板定义的 STEP_LIST。
            msg = f'Input step not found: {value}'
            msg += f'\nAvailable steps: {step_list}'
            raise ValueError(msg)

    # --dostep 优先级最高：忽略 --start/--end，只执行指定步骤。
    # check: ignore --start/end input if --dostep is specified
    if inps.doStep:
        inps.startStep = inps.doStep
        inps.endStep = inps.doStep

    # get list of steps to run
    # list.index(x) 返回 x 在列表中的位置编号，用它确定起始步骤和结束步骤的位置。
    idx0 = step_list.index(inps.startStep)
    idx1 = step_list.index(inps.endStep)
    if idx0 > idx1:
        # 起始步骤不能位于结束步骤之后。
        msg = f'start step "{inps.startStep}" is CAN NOT be after the end step "{inps.endStep}"'
        raise ValueError(msg)
    # Python 切片 step_list[idx0:idx1+1] 表示从 idx0 到 idx1，包括 idx1。
    run_steps = step_list[idx0:idx1+1]

    # 生成模板或仅绘图时，不执行处理步骤。
    # empty the step list if:
    # a) -g OR
    # b) iargs == ['--plot']
    if inps.generate_template or plot_only:
        run_steps = []

    # print mssage - processing steps
    if len(run_steps) > 0:
        # 多步骤运行时打印 logo，单步骤运行时打印紧凑版本信息。
        # for single step - compact version info
        if len(run_steps) == 1:
            print(mintpy.version.version_description)
        else:
            print(mintpy.version.logo)

        # datetime.datetime.now() 返回当前日期时间；os.getcwd() 返回当前工作目录。
        print(f'--RUN-at-{datetime.datetime.now()}--')
        print(f'Current directory: {os.getcwd()}')
        # os.path.basename(__file__) 只取当前脚本文件名，用于日志输出。
        print(f'Run routine processing with {os.path.basename(__file__)} on steps: {run_steps}')
        print(f'Remaining steps: {step_list[idx0+1:]}')
    print('-'*50)

    # 返回最终步骤列表，主工作流会逐个执行这些步骤。
    return run_steps


##########################################################################
def main(iargs=None):
    # CLI 入口：解析参数后调用主工作流实现。
    # parse
    inps = cmd_line_parse(iargs)

    # 延迟导入主工作流，减少只看帮助/版本时的启动开销。
    # import
    # from ... import ... 表示从另一个 Python 文件/模块里拿到一个函数。
    # 这里导入的是 src/mintpy/smallbaselineApp.py 里的 run_smallbaselineApp()。
    from mintpy.smallbaselineApp import run_smallbaselineApp

    # 将命令行参数对象交给 src/mintpy/smallbaselineApp.py 中的主流程执行。
    # run
    # 这里是 CLI 文件和主工作流文件的连接点：CLI 只负责解析参数，真正处理数据在 run_smallbaselineApp() 中。
    run_smallbaselineApp(inps)


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
