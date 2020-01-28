from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import StringField, IntegerField
from wtforms.validators import DataRequired
from wtforms.validators import ValidationError

def my_length_check(form, field):
    if len(field.data) < 2:
        raise ValidationError('Search string must be larger than 1 character')

def is_alpha_numeric(form, field):
    if not field.data.isalnum():
        raise ValidationError('Non alphanumeric characters are not allowed')


class Search_string(FlaskForm):
    search_string = StringField('', validators=[DataRequired(), my_length_check, is_alpha_numeric ])


class Inhibitor_used(FlaskForm):
    inhibitor = StringField('Used Inhibitor', validators=[DataRequired()])


class Threshold_pval(FlaskForm):
     threshold_pval = IntegerField('P-value Threshold')

class Threshold_foldchange(FlaskForm):
     threshold_foldchange = IntegerField('Foldchange Threshold')


class UploadForm(FlaskForm):
    uploaded_file = FileField('', validators=[FileAllowed(['tsv'], 'Should be ".tsv"')])
