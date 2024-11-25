import os
import stat
import glob


def add_project_root_to_bashrc():
    # 获取当前脚本的目录作为项目根目录
    project_root = os.path.abspath(os.path.dirname(__file__))
    bashrc_path = os.path.expanduser("~/.bashrc")

    # 在 .bashrc 中添加 PATH 设置
    path_export = f'\n# Add HiTE to PATH\nexport PATH={project_root}/module:{project_root}/tools:$PATH\n'

    # 检查 .bashrc 是否已经包含此路径
    with open(bashrc_path, 'r') as file:
        bashrc_content = file.read()

    if path_export.strip() not in bashrc_content:
        with open(bashrc_path, 'a') as file:
            file.write(path_export)
        print(f"Added {project_root} to PATH in ~/.bashrc.")
    else:
        print(f"{project_root} is already in PATH in ~/.bashrc.")

    # 提示用户手动加载 .bashrc 或打开新的终端
    print("To apply changes, run: source ~/.bashrc or open a new terminal.")


def add_execute_permission(file_paths):
    for file_path in file_paths:
        if os.path.isfile(file_path):  # 检查是否为文件
            st = os.stat(file_path)
            os.chmod(file_path, st.st_mode | stat.S_IEXEC)  # 添加执行权限
            print(f"Added execute permission to: {file_path}")


def main():
    # 获取当前脚本所在目录作为项目根目录
    project_root = os.path.abspath(os.path.dirname(__file__))

    # 将项目根目录添加到 PATH 并更新 .bashrc
    add_project_root_to_bashrc()

    # 定义需要添加执行权限的文件路径
    file_patterns = [
        os.path.join(project_root, "tools/*"),
        os.path.join(project_root, "module/*"),
        os.path.join(project_root, "bin/NeuralTE/tools/*"),
        os.path.join(project_root, "bin/LTR_FINDER_parallel-master/bin/LTR_FINDER.x86_64-1.0.7/ltr_finder"),
        os.path.join(project_root, "bin/HybridLTR-main/tools/*")
    ]

    # 为每个文件添加执行权限
    for pattern in file_patterns:
        for file in glob.glob(pattern):  # 使用 glob 处理通配符
            add_execute_permission([file])

    # 对LtrDetector进行重新编译
    print(f"Recompile LtrDetector...")
    project_root = os.path.abspath(os.path.dirname(__file__))
    LtrDetector_dir = project_root + '/bin/HybridLTR-main/bin/LtrDetector'
    LtrDetector_exe = LtrDetector_dir + '/bin/LtrDetector'
    if os.path.exists(LtrDetector_dir):
        recompile_cmd = 'cd ' + LtrDetector_dir + ' && rm -rf bin && cd src && make bin && make tr -j'
        os.system(recompile_cmd)
        if os.path.exists(LtrDetector_exe):
            print(f"finish recompile LtrDetector: {LtrDetector_exe}")
        else:
            print(f"fail to recompile LtrDetector: {LtrDetector_exe}")


if __name__ == "__main__":
    main()
