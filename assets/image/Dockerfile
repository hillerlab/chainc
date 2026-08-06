# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm AS builder

WORKDIR /app

COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release --all-features --bin chainc --locked && \
    strip target/release/chainc

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/chainc /usr/local/bin/chainc

# Set up non-root user
RUN useradd -m -u 1000 cuser && \
    chmod +x /usr/local/bin/chainc

USER cuser
WORKDIR /data

RUN chainc --help
