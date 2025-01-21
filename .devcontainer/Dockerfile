FROM mcr.microsoft.com/devcontainers/rust:1-1-bookworm

# Install development tools
RUN apt-get update \
  && export DEBIAN_FRONTEND=noninteractive \
  && apt-get -y install --no-install-recommends \
    python3-pip \
    python3-venv \
    python3-jinja2 \
    python3-pygments \
    perl \
    libssl-dev \
    libopenmpi-dev \
    libtbb-dev \
    libclang-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN \
  # sshd
  export DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  && apt-get -y install --no-install-recommends \
    openssh-server \
  # sshd_config
  && printf '%s\n' \
    'PermitRootLogin yes' \
    'PasswordAuthentication yes' \
    'PermitEmptyPasswords yes' \
    'UsePAM no' \
    > /etc/ssh/sshd_config.d/auth.conf \
  # ssh_config
  && printf '%s\n' \
    'Host *' \
    '    StrictHostKeyChecking no' \
    > /etc/ssh/ssh_config.d/ignore-host-key.conf \
  # delete passwd
  && passwd -d vscode \
  # Clean up
  && apt-get autoremove -y \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# RUN \
#   # Common packages
#   export DEBIAN_FRONTEND=noninteractive \
#   && apt-get update \
#   && apt-get -y install --no-install-recommends \
#   # Clean up
#   && apt-get autoremove -y \
#   && apt-get clean -y \
#   && rm -rf /var/lib/apt/lists/*
