from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed, FileRequired
from wtforms import StringField, IntegerField, DecimalField, SelectField
from wtforms.validators import DataRequired
from wtforms.validators import ValidationError

# valid if length of input string > 2 
def my_length_check(form, field):
    if len(field.data) < 2:
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
    chromosome = SelectField('Chromosome', \
                             choices = [('1,NC_000001.11','1'), ('2,NC_000002.12','2'), ('3,NC_000003.12','3'),\
                                        ('4,NC_000004.12','4'), ('5,NC_000005.10','5'), ('6,NC_000006.12','6'),\
                                        ('7,NC_000007.14','7'), ('8,NC_000008.11','8'), ('9,NC_000009.12','9'), \
                                        ('10,NC_000010.11','10'), ('11,NC_000011.10','11'), ('12,NC_000012.12','12'), \
                                        ('13,NC_000013.11','13'), ('14,NC_000014.9','14'), ('15,NC_000015.10','15'),\
                                        ('16,NC_000016.10','16'), ('17,NC_000017.11','17'), ('18,NC_000018.10','18'),\
                                        ('19,NC_000019.10','19'), ('20,NC_000020.11','20'), ('21,NC_000021.9','21'),\
                                        ('22,NC_000022.11','22'), ('X,NC_000023.11','X'), ('Y,NC_000024.10','Y')])
    genomic_loc_start = IntegerField('Genomic location Start', validators = [DataRequired(),is_positive])
    genomic_loc_end = IntegerField('Genomic Location End', validators = [DataRequired(),is_positive])

