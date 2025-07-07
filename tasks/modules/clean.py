from invoke import task

from ..config import (
    CLEAN_EXPRESSIONS,
)

### User config examples

# SNIPPET: expecting something like this
# CLEAN_EXPRESSIONS = [
#     "\"*~\"",
# ]

@task
def ls(cx):

    for clean_expr in CLEAN_EXPRESSIONS:
        cx.run(f'find . -type f -name {clean_expr} -print')

@task(pre=[ls], default=True)
def clean(cx):

    print("Deleting Targets")
    for clean_expr in CLEAN_EXPRESSIONS:
        cx.run(f'find . -type f -name {clean_expr} -delete')

