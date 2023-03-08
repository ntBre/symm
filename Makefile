TESTFLAGS = --test-threads=1 --nocapture
ARGS =

clippy:
	cargo clippy --all

test:
	cargo test -- ${TESTFLAGS} ${ARGS}
