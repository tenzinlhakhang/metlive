"""empty message

Revision ID: 210f39df5005
Revises: e7dcbdc2dd68
Create Date: 2019-10-28 22:04:47.700969

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '210f39df5005'
down_revision = 'e7dcbdc2dd68'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('user',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('username', sa.String(length=64), nullable=True),
    sa.Column('password_hash', sa.String(length=128), nullable=True),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_user_username'), 'user', ['username'], unique=True)
    op.drop_table('news_stories')
    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('news_stories',
    sa.Column('id', sa.BIGINT(), autoincrement=False, nullable=True),
    sa.Column('title', sa.VARCHAR(length=128), autoincrement=False, nullable=True),
    sa.Column('summary', sa.VARCHAR(length=256), autoincrement=False, nullable=True),
    sa.Column('story', sa.TEXT(), autoincrement=False, nullable=True)
    )
    op.drop_index(op.f('ix_user_username'), table_name='user')
    op.drop_table('user')
    # ### end Alembic commands ###
