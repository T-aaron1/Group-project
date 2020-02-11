from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import StringField, IntegerField, DecimalField
from wtforms.validators import DataRequired
from wtforms.validators import ValidationError

# valid if length of input string > 2 
def my_length_check(form, field):
    if len(field.data) < 2:
        raise ValidationError('Search string must be larger than 1 character')

# valid if length of input string equals 1 , used for chromosome
def length_one(form, field):
    if len(field.data) == 1:
        raise ValidationError('Search string must be larger than 1 character')
    
# valid if string contains only alphanumeric values
def is_alpha_numeric(form, field):
    if not field.data.isalnum():
        raise ValidationError('Non alphanumeric characters are not allowed')

# valid if input is positive
def is_positive(form, field):
    if field.data <= 0:
        raise ValidationError('0 or negative values are not allowed')

# FORM: used for the kinases    
class Search_string(FlaskForm):
    search_string = StringField('', validators=[DataRequired(), my_length_check, is_alpha_numeric ])

# FORM: used for the inhibitors
class Inhibitors(FlaskForm):
    inhibitor_name = StringField('', validators=[DataRequired(), my_length_check ])
    
# FORM: used for the uploaded file for the
#       phosphoproteomics analysis
class UploadForm(FlaskForm):
    threshold_pval = DecimalField('Log10 P-value Threshold', validators = [DataRequired(),is_positive])
    threshold_foldchange = DecimalField('Log2 Foldchange Threshold', validators = [DataRequired(),is_positive])
    cv_treatment_threshold = DecimalField('CV treatment Threshold', validators = [DataRequired(),is_positive])
    uploaded_file = FileField('', validators=[FileAllowed(['tsv'], 'Should be ".tsv"')])

# FORM: used for the chromosome coordinate search for phosphosites
class Phosphosite(FlaskForm):
    chromosome = StringField('Chromosome', validators=[DataRequired(), my_length_check, is_alpha_numeric,length_one ])
    genomic_loc_start = IntegerField('Genomic location Start')
    genomic_loc_end = IntegerField('Genomic Location End')

