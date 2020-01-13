from wtforms import Form, validators
from wtforms import FileField, StringField

class TmpForm(Form):
    uploaded_file = FileField('Archivo', [validators.Required(message = "Necesario"),
                                              validators.Regexp('.+\.csv$', message = "no es csv")
                                              
    ]
    )


class UploadPhosphoproteomics(Form):
    uploaded_file = FileField('Phosphoproteomics',[ validators.Required(message = "No file selected"),
                                                      validators.Regexp('.+\.tsv$', message = 'Not a ".tsv" file')]
    )
