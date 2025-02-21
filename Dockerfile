FROM nvidia/cuda:11.1.1-devel-ubuntu20.04
MAINTAINER Comma "douhao_zy@163.com"

# Set an encoding to make things work smoothly.
ENV LANGUAGE en_US.UTF-8
ENV LC_ALL C
ENV LANG en_US.UTF-8
ENV TZ US/Pacific
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility
ENV PATH /usr/local/cuda-11.1/bin:$PATH
ENV LD_LIBRARY_PATH /usr/local/cuda-11.1/lib64:$LD_LIBRARY_PATH
ENV CUDA_HOME /usr/local/cuda-11.1
ARG DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]

# 安装依赖
RUN sed -i 's|http://archive.ubuntu.com/ubuntu/|https://mirrors.ustc.edu.cn/ubuntu/|g' /etc/apt/sources.list \
 && set -ex \
 && apt-get update \
 && apt-get install -y --fix-missing \
    g++ \
    cmake \
    gfortran \
    unzip \
    curl \
    libgl1 \
    libgdal-dev \
    libgomp1 \
    gnuplot \
    gmt \
    gmt-dcw \
    gmt-gshhg \
    ghostscript \
    gdal-bin \
    pkg-config \
    libfftw3-dev \
    libmotif-dev \
    libopencv-dev \
    libx11-dev \
    libeccodes-dev \
 && echo done

# 编译安装MintPy
COPY ./Miniconda3-latest-Linux-x86_64.sh /opt/Miniconda3-latest-Linux-x86_64.sh
COPY ./MintPy /opt/MintPy
RUN bash /opt/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/Miniconda3 \
 && eval "$(/opt/Miniconda3/bin/conda shell.bash hook)" \
 && conda init \
 && conda update -n base conda \
 && conda update -y --all \
 && conda install -c conda-forge mamba \
 && mamba install -c conda-forge -y --file /opt/MintPy/requirements.txt \
 && chmod 777 -R /opt/MintPy/* \
 && sed -i '6s/^/#/g' /root/.bashrc \
 && echo done

# 编译安装ISCE-2.6.2
COPY ./isce2-2.6.2 /opt/isce2
RUN source /root/.bashrc \
 && cd /opt/isce2 \
 && mkdir build \
 && cd build \
 && cmake .. -DCMAKE_INSTALL_PREFIX=/opt/isce/install/location \ 
 && make install -j4 \
 && echo done
COPY ./isce2-2.6.2/license.py /opt/isce/install/location/packages/isce/license.py

# 编译安装MSBAS
COPY ./MSBAS-6 /opt/MSBAS
RUN source /root/.bashrc \
 && cd /opt/MSBAS/msbasv6 && make -j4 \
 && cd /opt/MSBAS/msbas_extract && make -j4 \
 && echo done
 
# 配置文件拷贝
COPY ./config_files/bgyr.cpt /usr/share/gmt/cpt/bgyr.cpt
COPY ./config_files/netrc /root/.netrc
COPY ./config_files/cdsapirc /root/.cdsapirc
COPY ./config_files/dask.yaml /root/.config/dask/dask.yaml
COPY ./config_files/smallbaselineSentinel.txt /root/smallbaselineSentinel.txt

# 环境设置
ENV PYTHONPATH /opt/isce/install/location/packages
ENV ISCE_HOME /opt/isce/install/location/packages/isce
ENV PATH $ISCE_HOME/applications:$PATH
ENV PATH /opt/isce2/contrib/stack/topsStack:$PATH
ENV MINTPY_HOME=/opt/MintPy
ENV PATH $MINTPY_HOME/src/mintpy/cli:$PATH
ENV PYTHONPATH ${PYTHONPATH}:$MINTPY_HOME/src:$PATH
ENV VRT_SHARED_SOURCE 0
ENV HDF5_DISABLE_VERSION_CHECK 2
ENV HDF5_USE_FILE_LOCKING FALSE
ENV MSBAS_HOME /opt/MSBAS
ENV PATH $MSBAS_HOME/msbasv6:$PATH
ENV PATH $MSBAS_HOME/msbas_extract:$PATH

# 解决一个库的版本问题
CMD rm /usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
 && ln -s /opt/Miniconda3/lib/libstdc++.so.6.0.33 /usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
 && /bin/bash
