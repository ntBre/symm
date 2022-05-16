TESTFLAGS = --test-threads=1 --nocapture

test:
	cargo test -- ${TESTFLAGS}
