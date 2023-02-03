from docutils import nodes
from docutils.nodes import Element, Node
from sphinx.util.docutils import SphinxDirective
from docutils.parsers.rst import directives
from docutils.parsers.rst.directives.admonitions import BaseAdmonition
from typing import Any, Callable, Dict, List
from sphinx.writers.html import HTMLTranslator
OptionSpec = Dict[str, Callable[[str], Any]]
from sphinx.locale import _

class exercise_node(nodes.Admonition, nodes.Element):
    pass

class SSHExercise(BaseAdmonition, SphinxDirective):
    node_class = exercise_node
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec: OptionSpec = {
        'class': directives.class_option,
        'name': directives.unchanged,
    }

    def run(self) -> List[Node]:
        if not self.options.get('class'):
            self.options['class'] = ['admonition-important']

        (todo,) = super().run()  # type: Tuple[Node]
        if isinstance(todo, nodes.system_message):
            return [todo]
        elif isinstance(todo, exercise_node):
            todo.insert(0, nodes.title(text=_('Exercise')))
            todo['docname'] = self.env.docname
            self.add_name(todo)
            self.set_source_info(todo)
            self.state.document.note_explicit_target(todo)
            return [todo]
        else:
            raise RuntimeError  # never reached here

    # def run(self):
    #     paragraph_node = nodes.paragraph(text='Hello World')
    #     return [paragraph_node]


def visit_exercise_node(self: HTMLTranslator, node: exercise_node) -> None:
    if self.config.todo_include_todos:
        self.visit_admonition(node)
    else:
        raise nodes.SkipNode


def depart_exercise_node(self: HTMLTranslator, node: exercise_node) -> None:
    self.depart_admonition(node)

class solution_node(nodes.Admonition, nodes.Element):
    pass

class SSHSolution(BaseAdmonition, SphinxDirective):
    node_class = solution_node
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec: OptionSpec = {
        'class': directives.class_option,
        'name': directives.unchanged,
    }

    def run(self) -> List[Node]:
        if not self.options.get('class'):
            self.options['class'] = ['admonition-important']

        (todo,) = super().run()  # type: Tuple[Node]
        if isinstance(todo, nodes.system_message):
            return [todo]
        elif isinstance(todo, solution_node):
            todo.insert(0, nodes.title(text=_('Solution')))
            todo['docname'] = self.env.docname
            self.add_name(todo)
            self.set_source_info(todo)
            self.state.document.note_explicit_target(todo)
            return [todo]
        else:
            raise RuntimeError  # never reached here

    # def run(self):
    #     paragraph_node = nodes.paragraph(text='Hello World')
    #     return [paragraph_node]


def visit_solution_node(self: HTMLTranslator, node: solution_node) -> None:
    if self.config.todo_include_todos:
        self.visit_admonition(node)
    else:
        raise nodes.SkipNode


def depart_solution_node(self: HTMLTranslator, node: solution_node) -> None:
    self.depart_admonition(node)

class prereq_node(nodes.Admonition, nodes.Element):
    pass

class Prereq(BaseAdmonition, SphinxDirective):
    node_class = prereq_node
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec: OptionSpec = {
        'class': directives.class_option,
        'name': directives.unchanged,
    }

    def run(self) -> List[Node]:
        if not self.options.get('class'):
            self.options['class'] = ['admonition-important']

        (todo,) = super().run()  # type: Tuple[Node]
        if isinstance(todo, nodes.system_message):
            return [todo]
        elif isinstance(todo, prereq_node):
            todo.insert(0, nodes.title(text=_('Prerequesites')))
            todo['docname'] = self.env.docname
            self.add_name(todo)
            self.set_source_info(todo)
            self.state.document.note_explicit_target(todo)
            return [todo]
        else:
            raise RuntimeError  # never reached here

    # def run(self):
    #     paragraph_node = nodes.paragraph(text='Hello World')
    #     return [paragraph_node]


def visit_prereq_node(self: HTMLTranslator, node: prereq_node) -> None:
    if self.config.todo_include_todos:
        self.visit_admonition(node)
    else:
        raise nodes.SkipNode


def depart_prereq_node(self: HTMLTranslator, node: prereq_node) -> None:
    self.depart_admonition(node)

def setup(app):

    app.add_node(solution_node, html=(visit_solution_node, depart_solution_node))
    app.add_node(exercise_node, html=(visit_exercise_node, depart_exercise_node))
    app.add_node(prereq_node, html=(visit_prereq_node, depart_prereq_node))

    app.add_directive('sshsolution', SSHSolution)
    app.add_directive('sshexercise', SSHExercise)
    app.add_directive('prereq', Prereq)

    return {
    'version': '0.1',
    'parallel_read_safe': True,
    'parallel_write_safe': True,
    }
