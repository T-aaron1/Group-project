from wtforms import Form, validators
from wtforms import FileField, StringField


class UploadPhosphoproteomics(Form):
    uploaded_file = StringField('Phosphoproteomics',[
                               validators.length(min=4, max=25, message = 'miaosd')]
    )
