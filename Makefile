.PHONY: verify negative-control reproduce

verify:
	python3 code/verify_sis.py
	python3 -O code/verify_sis.py

negative-control:
	python3 code/test_negative_control.py

reproduce: verify negative-control
