#!/bin/bash

# FTP服务器信息
ftp_server="ftp://human.big.ac.cn/"
username="leishi@lzu.edu.cn"
password="Daxia12345!"

# 目标文件夹
target_folder="HRA001130"

# 待下载的文件夹列表
folders=("HRR338285" "HRR338284" "HRR338086" "HRR338087")

# 创建下载目录
download_dir="downloads"
mkdir -p "$download_dir"

# 登录FTP服务器
wget --user="$username" --password="$password" --recursive --no-clobber --no-parent --no-directories --level=1 --directory-prefix="$download_dir" "$ftp_server$target_folder/"

# 多线程下载各个文件夹
for folder in "${folders[@]}"; do
  {
    wget --user="$username" --password="$password" --recursive --no-clobber --no-parent --no-directories --level=1 --directory-prefix="$download_dir/$folder" "$ftp_server$target_folder/$folder/"
  } &
done

# 等待所有后台进程完成
wait

echo "下载完成"
