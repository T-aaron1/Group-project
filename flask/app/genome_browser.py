from flask import Flask, Blueprint, render_template, redirect, request, url_for, session, flash
import forms
from flask_wtf import CsrfProtect
import os
import pandas as pd
import queries




genome_browser_blueprint = Blueprint(
    'genome_browser', __name__,
    template_folder = 'genome_browser'
)



@genome_browser_blueprint.route('/<chromosome>')
def genome_viewer_chrom(chromosome):
    pass
