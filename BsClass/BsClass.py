import os
import sys
import argparse
from snakemake import snakemake

class BsClass:
    def __init__(self, workflowfile, cores, configfile):
        self.workflowfile = workflowfile
        self.cores = cores
        self.configfile = configfile

    def run_workflow(self):
        # 使用 Snakemake 运行工作流
        try:
            success = snakemake(
                snakefile=self.workflowfile,  # 使用提供的工作流文件
                cores=self.cores,  # 根输入调整核心数
                config={"default_configfile": self.configfile},  # 将YAML配置文件路径传递给Snakemake
                keepgoing=True,  # 对应 --keep-going
                printshellcmds=True,  # 对应 --printshellcmds
                reason=True,  # 对应 --reason
                use_conda=True,  # 对应 --use-conda
                rerun_incomplete=True  # 对应 --rerun-incomplete
            )
            if not success:
                raise Exception("Snakemake workflow failed.")
        except Exception as e:
            print(f"Error running workflow: {e}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="BsClass: A tool for NGS data analysis")
    
    # 允许选择 workflow (NGS, Pacbio, ONT)
    parser.add_argument('--workflow', required=True, choices=['NGS', 'Pacbio', 'ONT'], help="Specify the workflow type (NGS, Pacbio, ONT)")
    
    # 新增选择核心数的参数
    parser.add_argument('--cores', type=int, default=4, help="Number of CPU cores to use (default: 4)")
    
    # 新增 YAML 配置文件路径参数
    parser.add_argument('--configfile', required=True, help="Path to the YAML configuration file")
    
    args = parser.parse_args()

    # 获取当前 Python 脚本所在的目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 将选择的 workflow 映射到相应的文件路径
    workflow_map = {
        'NGS': os.path.join(script_dir, 'workflow', 'NGS.smk'),
        'Pacbio': os.path.join(script_dir, 'workflow', 'Pacbio.smk'),
        'ONT': os.path.join(script_dir, 'workflow', 'ONT.smk')
    }
    
    # 确定工作流文件路径
    workflowfile = workflow_map[args.workflow]

    # 创建BsClass实例并运行工作流
    bs = BsClass(workflowfile=workflowfile, cores=args.cores, configfile=args.configfile)
    bs.run_workflow()

if __name__ == "__main__":
    main()

