from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import StringField, IntegerField, DecimalField
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


class UploadForm(FlaskForm):
    threshold_pval = DecimalField('Log10 P-value Threshold')
    threshold_foldchange = DecimalField('Log2 Foldchange Threshold')
    uploaded_file = FileField('', validators=[FileAllowed(['tsv'], 'Should be ".tsv"')])

class Phosphosite(FlaskForm):
    chromosome = StringField('Chromosome', validators=[DataRequired(), my_length_check, is_alpha_numeric ])
    genomic_loc_start = IntegerField('Genomic location Start')
    genomic_loc_end = IntegerField('Genomic Location End')

class Inhibitors(FlaskForm):
    inhibitor_name = StringField('', validators=[DataRequired(), my_length_check ])


