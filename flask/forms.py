from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired




class UploadForm(FlaskForm):
    uploaded_file = FileField('Phosphoproteomics', validators=[
        FileAllowed(['tsv'], 'Should be ".tsv"')
    ])
