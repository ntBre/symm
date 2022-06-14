TESTFLAGS = --test-threads=1 --nocapture
ARGS =

test:
	cargo test -- ${TESTFLAGS} ${ARGS}
