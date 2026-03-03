SHELL := /bin/bash

.PHONY: check py sync run

NOTEBOOK_FIND = find . \
	-mindepth 1 -type d \( -name '.*' -o -name '.ipynb_checkpoints' \) -prune -o \
	-type f -name '*.ipynb' -print0

check:
	@paired=0; unpaired=0; total=0; \
	while IFS= read -r -d '' nb; do \
		total=$$((total + 1)); \
		py_file="$${nb%.ipynb}.py"; \
		if [[ -f "$$py_file" ]]; then \
			echo "PAIRED   $$nb <-> $$py_file"; \
			paired=$$((paired + 1)); \
		else \
			echo "UNPAIRED $$nb"; \
			unpaired=$$((unpaired + 1)); \
		fi; \
	done < <(find . -mindepth 1 -type d \( -name '.*' -o -name '.ipynb_checkpoints' \) -prune -o -type f -name '*.ipynb' -print0); \
	echo; \
	echo "Summary: total=$$total paired=$$paired unpaired=$$unpaired"

py:
	@$(NOTEBOOK_FIND) | while IFS= read -r -d '' nb; do \
		py_file="$${nb%.ipynb}.py"; \
		if [[ ! -f "$$py_file" ]]; then \
			echo "Pairing $$nb -> $$py_file"; \
			jupytext --set-formats ipynb,py:percent "$$nb"; \
			jupytext --sync "$$nb"; \
		fi; \
	done

sync:
	@$(NOTEBOOK_FIND) | while IFS= read -r -d '' nb; do \
		echo "Syncing $$nb"; \
		jupytext --sync "$$nb"; \
	done

run:
	@$(NOTEBOOK_FIND) | while IFS= read -r -d '' nb; do \
		echo "Executing $$nb"; \
		jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 "$$nb"; \
	done
