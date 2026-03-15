FROM ubuntu:22.04

# Install build tools
RUN apt update && apt install -y \
    build-essential \
    cmake \
    git \
    libboost-all-dev

# Create work directory
WORKDIR /enumhyp

# Clone enumhyp
RUN git clone https://github.com/goodefroi/enumhyp.git /enumhyp

# Build enumhyp
RUN mkdir -p build && cd build && cmake .. && make -j4


# Set default app (Linux output goes to /enumhyp/build/bin/enumhyp)
ENTRYPOINT ["/enumhyp/build/bin/enumhyp"]

