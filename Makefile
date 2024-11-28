run_test_circ:
	python -m unittest tests/UpUCCSDG_test.py

prod:
	python -m unittest tests/paulistring_test.py

my_circ:
	python -m unittest tests/qc_tests.py
