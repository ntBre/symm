TESTFLAGS = --test-threads=1 --nocapture
ARGS =

clippy:
	cargo clippy --tests

test:
	cargo test -- ${TESTFLAGS} ${ARGS}
