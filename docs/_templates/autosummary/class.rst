{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree: .
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: .
   {% for item in methods %}
      {%- if item != '__init__' %}
      {{ objname }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}