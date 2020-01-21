from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import StringField
from wtforms.validators import DataRequired

from wtforms.fields import IntegerField


class Inhibitor_used(FlaskForm):
    inhibitor = StringField('Used Inhibitor', validators=[DataRequired()])

class Threshold_pval(FlaskForm):
     threshold_pval = IntegerField('P-value Threshold', validators=[DataRequired()])

class Threshold_foldchange(FlaskForm):
     threshold_foldchange = IntegerField('Foldchange Threshold', validators=[DataRequired()])


class UploadForm(FlaskForm):
    uploaded_file = FileField('', validators=[
        FileAllowed(['tsv'], 'Should be ".tsv"')
    ])
