{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions

   .. autosummary::
      :toctree: generated/
      {% for item in functions %}
         {{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes

   .. autosummary::
      :toctree: generated/
      {% for item in classes %}
      {% if item[-4:] == "Test" %}
      :template: autosummary/inherits_TestCase_class.rst
      {% else %}
      :template: autosummary/class.rst
      {% endif %}
         {{ item }}

      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions

   .. autosummary::
      :toctree: generated/
      {% for item in exceptions %}
         {{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}
